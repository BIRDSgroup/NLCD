import os 
os.environ['OMP_NUM_THREADS'] = '1' 
from tensorflow_probability import distributions as tfd
from tensorflow.keras.layers import Dense
import tensorflow as tf
import matplotlib.pyplot as plt
from scipy.stats import norm
from tensorflow.keras.models import Sequential
import random
import numpy as np
from sklearn.svm import SVR
import pandas as pd
from sklearn.kernel_ridge import KernelRidge
from multiprocessing import Pool
from numpy.random import SeedSequence
from sklearn.metrics import mean_squared_error
import sys

def gnll_eval(y,alpha, mu, sigma):
    """ Computes the mean negative log-likelihood loss of y given the mixture parameters.
    """
    gm = tfd.MixtureSameFamily(
        mixture_distribution=tfd.Categorical(probs=alpha),
        components_distribution=tfd.Normal(
            loc=mu,       
            scale=sigma))
    log_likelihood = gm.log_prob(tf.transpose(y))
    return -tf.reduce_mean(log_likelihood, axis=-1)


def calculate_pvalue(original,loss_list,greater=False):
    '''
    calculate the p value 
    '''
    if(greater==False):
        return sum(i < original for i in loss_list)/len(loss_list)
    else:
        #It is reverse in test 2
        return sum(i >= original for i in loss_list)/len(loss_list)


def get_prob(L,A):
    '''
    Function to calculate the overlap score 
    '''
    norm_const=0
    probs=[]
    unique_values = np.unique(L)
    for value in unique_values:
        indices = np.where(L == value)[0]
        mu = np.mean(A[indices])
        sigma = np.std(A[indices])
        p_L = np.count_nonzero(L == value) / len(L)
        p_AgivenL = norm.pdf(A, mu, sigma)
        p_LgivenA = p_AgivenL * p_L
        probs.append(p_LgivenA)
        norm_const += p_LgivenA

    diff=np.minimum.reduce(probs)


    return diff/norm_const


def FI_score(x,y,overlap ):
    '''
    Function to calcuate the overlap score 
    '''

    # Iterate over the pairs of lists and calculate element-wise absolute differences
    diff = np.abs(np.subtract(x, y))

    return np.sum(np.multiply(diff,overlap))/len(x)


def compute_4_loss(L,A,B,algo):
    '''
    Function to calculate the predictions of y for different values of L in the 4th test 
    '''
    unique_values = np.unique(L)
    y_pred = []

    if algo == "SVR":
        regressor = SVR(kernel='rbf')
        X = np.column_stack((L, A))
        regressor.fit(X, B)

        for value in unique_values:
            L_1 = np.full_like(L, value)
            X_ = np.column_stack((L_1, A))
            y_predict_ = regressor.predict(X_)
            y_pred.append(y_predict_)

    
    elif algo=="ANN":
        model = Sequential()
        model.add(Dense(12, input_shape=(2,), activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='linear'))

        X = np.column_stack((L, A))

        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_squared_error'])
        model.fit(X, B, epochs=150, batch_size=64,verbose=0)

        for value in unique_values:
            L_1 = np.full_like(L, value)
            X_ = np.column_stack((L_1, A))
            y_predict_ = model.predict(X_).reshape(-1)
            y_pred.append(y_predict_)

    elif algo == "KRR":
        regressor = KernelRidge(kernel='rbf')
        X = np.column_stack((L, A))
        regressor.fit(X, B)

        for value in unique_values:
            L_1 = np.full_like(L, value)
            X_ = np.column_stack((L_1, A))
            y_predict_ = regressor.predict(X_)
            y_pred.append(y_predict_)
    else:
        print("Invalid Algorithm")

    return y_pred


def stratify_permute_variable(L, variable):
    '''
    Function to stratify permute a variable 
    '''
    unique_values = np.unique(L)
    permuted_variable = np.empty_like(variable)

    for value in unique_values:
        indices = np.where(L == value)[0]
        permuted_indices = random.sample(list(indices), len(indices))
        permuted_variable[indices] = variable[permuted_indices]

    return permuted_variable


def test_4(L,A,B,shuffles,algo,test_2=False):
    '''
    Function for the 4th test , same function is used for the second test 
    '''

    overlap=get_prob(L,A)
    perm_loss=[]
    for i in range(shuffles):
            if test_2==False:
                y_pred=compute_4_loss(L,A,stratify_permute_variable(L,B),algo)
            else:
            
                y_pred=compute_4_loss(L,A,np.random.permutation(B),algo) 
            total_FI=0
            count=0
        
            for i in range(len(y_pred)):
                for j in range(i + 1, len(y_pred)):
                    total_FI+=FI_score(y_pred[i],y_pred[j],overlap)
                    count+=1
            total_FI/=count
            perm_loss.append(total_FI)
    

    y_pred_original=compute_4_loss(L,A,B,algo)
    original_loss=0
    count=0
    for i in range(len(y_pred_original)):
        for j in range(i + 1, len(y_pred_original)):
            original_loss+=FI_score(y_pred_original[i],y_pred_original[j],overlap)
            count+=1
    original_loss/=count
    return [calculate_pvalue(original_loss,perm_loss,test_2),np.sum(overlap)/len(overlap)]


def compute_1_loss(x_test,y_test):
    '''
    Function to calculate the loss in the first test using nll 
    '''
    unique_values = np.unique(x_test)

    y_mean = np.empty_like(y_test)
    y_std = np.empty_like(y_test)

    for value in unique_values:
        indices = np.where(x_test == value)[0]
        mu, sigma = np.mean(y_test[indices]), np.std(y_test[indices])
        y_mean[indices] = mu
        y_std[indices] = sigma
    #alpha is the probability of each distributiton, here we are all having it as 1 
    alpha=np.ones((len(y_mean),1))
    return gnll_eval(y_test,alpha,y_mean.reshape((-1,1)),y_std.reshape((-1,1))).numpy()


def test_1(L,B,shuffles):
    '''
    Function for the first test 
    '''
    perm_loss=[]
    for i in range(0,shuffles):
        perm_loss.append(compute_1_loss(L,np.random.permutation(B)))
    original_loss=compute_1_loss(L,B)
    return calculate_pvalue(original_loss,perm_loss)    


def test_3_loss(A,B,algo):
    '''
    function to return the loss in the third function 
    '''
    
    if algo == "SVR":
        regressor = SVR(kernel='rbf')
        regressor.fit(A.reshape(-1, 1), B)
        y_predict = regressor.predict(A.reshape(-1, 1))

    elif algo == "ANN":
        model = Sequential([
            Dense(12, input_shape=(1,), activation='relu'),
            Dense(8, activation='relu'),
            Dense(1, activation='linear')
        ])
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_squared_error'])
        model.fit(A, B, epochs=150, batch_size=64, verbose=0)
        y_predict = model.predict(A)

    elif algo == "KRR":
        regressor = KernelRidge(kernel='rbf')
        regressor.fit(A.reshape(-1, 1), B)
        y_predict = regressor.predict(A.reshape(-1, 1))

    mse = mean_squared_error(B, y_predict)

    return mse


def test_3(L,A,B,shuffles,algo):
    '''
    function for the third tes t
    '''

    unique_values = np.unique(L)

    p_values = []

    for value in unique_values:
        indices = np.where(L == value)
        A_value = A[indices]
        B_value = B[indices]
        
        original_loss = test_3_loss(A_value, B_value, algo)
        perm_loss = [test_3_loss(A_value, np.random.permutation(B_value), algo) for _ in range(shuffles)]
        p_value = calculate_pvalue(original_loss, perm_loss)
        p_values.append(p_value)

    max_p_value = max(p_values)

    return max_p_value


def test_2(L,A,B,shuffles,algo):         #using the regression method
    if(algo=="SVR"):  #SVM
        regressor = SVR(kernel = 'rbf')
        regressor.fit(B.reshape(-1,1),A)
        y_predict=regressor.predict(B.reshape(-1,1))
        y_resid=A-y_predict
    elif(algo=="KRR"):  #KRR
        regressor = KernelRidge(kernel = 'rbf')
        regressor.fit(B.reshape(-1,1),A)
        y_predict=regressor.predict(B.reshape(-1,1))
        y_resid=A-y_predict
    elif(algo=="ANN"):
        model = Sequential()
        model.add(Dense(12, input_shape=(1,), activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='linear'))
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_squared_error'])
        model.fit(B, A, epochs=150, batch_size=64,verbose=0)
        y_pred=model.predict(B)
        y_resid=A-y_pred.reshape((len(y_pred),))
    return test_1(L,y_resid,shuffles)

def test_2_resid(L,A,B,shuffles,algo):
    if(algo=="SVR"):  #SVM
        regressor = SVR(kernel = 'rbf')
        regressor.fit(B.reshape(-1,1),A)
        y_predict=regressor.predict(B.reshape(-1,1))
        y_resid=A-y_predict
    elif(algo=="KRR"):  #KRR
        regressor = KernelRidge(kernel = 'rbf')
        regressor.fit(B.reshape(-1,1),A)
        y_predict=regressor.predict(B.reshape(-1,1))
        y_resid=A-y_predict
    elif(algo=="ANN"):
        model = Sequential()
        model.add(Dense(12, input_shape=(1,), activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='linear'))
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_squared_error'])
        model.fit(B, A, epochs=150, batch_size=64,verbose=0)
        y_pred=model.predict(B)

        y_resid=A-y_pred.reshape((len(y_pred),))
    Astar=y_predict+np.random.permutation(y_resid)
    return test_4(L,Astar,B,shuffles,algo,True)

def pcorrection(pval,shuffles):
    return (pval*shuffles+1)/(shuffles+2)


def combine_tests(L,A,B,shuffles,algo):
    '''
    Function to combine all the tests 
    '''
    #LA_p=test_1(L,A,shuffles)
    LB_p=test_1(L,B,shuffles)
    LAgvnB,overlapscore1=test_4(L,B,A,shuffles,algo,True)
    ABgvnL=test_3(L,A,B,shuffles,algo)
    LindBgvnA,overlapscore2=test_4(L,A,B,shuffles,algo)
    # p value correction for all the tests 
    LB_p_corr=pcorrection(LB_p,shuffles)
    LAgvnB_corr=pcorrection(LAgvnB,shuffles)
    ABgvnL_corr=pcorrection(ABgvnL,shuffles)
    LindBgvnA_corr=pcorrection(LindBgvnA,shuffles)
    p_final=np.max([LB_p_corr,LAgvnB_corr,ABgvnL_corr,LindBgvnA])
    return [p_final,LB_p_corr,LAgvnB_corr,ABgvnL_corr,LindBgvnA_corr,overlapscore1,overlapscore2]


