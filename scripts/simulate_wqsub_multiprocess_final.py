
import time
#to keep track of the total time taken 
st = time.time()
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
import sys
# first argument is for the algorithm
algo=sys.argv[1]
# second argument is for the data file 
data= sys.argv[2] 
# third argument is for the output file 
outputname= sys.argv[3]
# fourth argument is for the number of permutations 
shuffles=int(sys.argv[4])


#read the input file 
fo=open(data, "r")
L=[]
A=[]
B=[]
# the simulated data is generated with a seed at the start of the file
# so we skip it, but not for the yeast data 
# if yeast is present then dont read the seed line
if(data.find('yeast')==-1): #if yeast is not present
    line=fo.readline() # read the line, the control will start from the next line 

while(1):
    # the first line of the trio is the parameter configuration, skip it 
    line=fo.readline()
    # read the genotype information 
    line=fo.readline()
    if line== "" :
        break
    l = [j for j in line.split()]
    L.append([int(i) for i in l])
    line=fo.readline()
    a = [j for j in line.split()]
    A.append([float(i) for i in a])
    line=fo.readline()
    b = [j for j in line.split()]
    B.append([float(i) for i in b])
#combine all the samples into a list 
dataset_linear = [i for i in zip(L,A,B)]
fo.close()

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



def compute_loss(x_test,y_test):
    """ Evaluate the model to get the loss by calculating the mle of mu and sigma 
    """
    #get the indices of the genotype values 0 and 1
    indices_1 = [i for i, x in enumerate(x_test) if x == 1]
    indices_0 = [i for i, x in enumerate(x_test) if x == 0]
    #calculate the mle of mu and sigma for 0 and 1
    mu_0=np.mean(y_test[indices_0])
    mu_1=np.mean(y_test[indices_1])
    sigma_0=np.std(y_test[indices_0])
    sigma_1=np.std(y_test[indices_1])
    #initialize the array and set the value 
    y_mean=np.zeros((len(y_test),1))
    y_mean[indices_1]=mu_1
    y_mean[indices_0]=mu_0
    y_std=np.zeros((len(y_test),1))
    y_std[indices_1]=sigma_1
    y_std[indices_0]=sigma_0
    #alpha is the probability of each distributiton, here we are all having it as 1  
    alpha=np.ones((len(y_mean),1))
    return gnll_eval(y_test,alpha,y_mean,y_std).numpy()
    
def reshapevar(X):
  """
  Function to reshape the vector for the input 
  """
  return np.reshape(X,newshape=(len(X),-1))

#used in L->A and L->B testing 
def shuffleBtimes(P,Q,B):
    """ Shuffle Q B times and compute the loss 
    """
    loss=[]
    for i in range(0,B):
      loss.append(compute_loss(P,np.random.permutation(Q)))
    return loss


def calculate_pvalue(original,loss_list):
    '''
    calculate the p value 
    '''
    return sum(i < original for i in loss_list)/len(loss_list)
    
    

def compute_third_testloss(L,A,B,algo):
    '''
    Function to calculate the Yprediction for L=0,A and L=1,A for different algorithms  
    '''
    if(algo=="SVR"):  
        regressor = KernelRidge(kernel = 'rbf')
        X=np.vstack((L,A)).T
        regressor.fit(X,B)
	#with L = 0 and L = 1 give the inputs 
        L_ones=np.ones((L.shape))
        L_minus=np.ones((L.shape))*0
        X_zero=np.vstack((L_minus,A)).T
        X_one=np.vstack((L_ones,A)).T
        y_predict_zero=regressor.predict(X_zero)
        y_predict_one=regressor.predict(X_one)
        return [y_predict_one,y_predict_zero]
    
    elif(algo=="ANN"):
        model = Sequential()
        model.add(Dense(12, input_shape=(2,), activation='relu'))
        model.add(Dense(8, activation='relu'))
        model.add(Dense(1, activation='linear'))
        X=np.vstack((L,A)).T
        model.compile(loss='mean_squared_error', optimizer='adam', metrics=['mean_squared_error'])
        model.fit(X, B, epochs=150, batch_size=64,verbose=0)
        L_ones=np.ones((L.shape))
        L_minus=np.ones((L.shape))*0
        X_zero=np.vstack((L_minus,A)).T
        X_one=np.vstack((L_ones,A)).T
        y_predict_zero=model.predict(X_zero)
        y_predict_one=model.predict(X_one)
        y_predict_zero=y_predict_zero.reshape((len(y_predict_zero,)))
        y_predict_one=y_predict_one.reshape((len(y_predict_one,)))
        return [y_predict_one,y_predict_zero]
    
    elif(algo=="KRR"):

        regressor = KernelRidge(kernel = 'rbf')
        X=np.vstack((L,A)).T
        regressor.fit(X,B)
        L_ones=np.ones((L.shape))
        L_minus=np.ones((L.shape))*0
        X_zero=np.vstack((L_minus,A)).T
        X_one=np.vstack((L_ones,A)).T
        y_predict_zero=regressor.predict(X_zero)
        y_predict_one=regressor.predict(X_one)
        return [y_predict_one,y_predict_zero]

    

def get_prob(L,A):
    '''
    this function will calculate the p(L|A) and returns the minimum of both 
    '''
    indices_1 = [i for i, x in enumerate(L) if x == 1]
    indices_0 = [i for i, x in enumerate(L) if x == 0]
    mu_0=np.mean(A[indices_0])
    mu_1=np.mean(A[indices_1])
    sigma_0=np.std(A[indices_0])
    sigma_1=np.std(A[indices_1])
    p_L0= np.count_nonzero(L==0)/len(L)
    p_L1= np.count_nonzero(L==1)/len(L)
    p_AgivenL0 = norm.pdf(A,mu_0,sigma_0)
    p_AgivenL1 = norm.pdf(A,mu_1,sigma_1)
    p_L0givenA = p_AgivenL0 * p_L0
    p_L1givenA = p_AgivenL1 * p_L1
    const= p_L0givenA+p_L1givenA
    p_L0givenA_norm= p_L0givenA/const
    p_L1givenA_norm=p_L1givenA/const
    diff=np.minimum(p_L0givenA_norm,p_L1givenA_norm)
    return diff

def stratify_B_n_times_diff(L,A,B,n,algo):
    '''
    this function computes the FI score for n permutations, returns a list of losses
    '''
    loss=[]
    indices_1 = [i for i, x in enumerate(L) if x == 1]
    indices_0 = [i for i, x in enumerate(L) if x == 0]
    diff=get_prob(L,A)
    for i in range(0,n):
        B_dist_temp=np.zeros(len(B))
        mod_indices_1=random.sample(indices_1,len(indices_1))
        for i in range(len(indices_1)):
            B_dist_temp[indices_1[i]]= B[mod_indices_1[i]]

        mod_indices_0=random.sample(indices_0,len(indices_0))
        for i in range(len(indices_0)):
            B_dist_temp[indices_0[i]]= B[mod_indices_0[i]]

        y_pred_ones,y_pred_zeros=compute_third_testloss(L,A,B_dist_temp,algo)
        
       
        loss.append(sum(np.multiply(abs(y_pred_zeros-y_pred_ones),diff))/len(y_pred_ones))
    
    return loss


def test_1(L,T,shuffles):
    '''
    this function conducts the first test, using the same function for the second test, returns a p value
    '''
    T_shuffle=np.copy(T)
    loss_list_LT=shuffleBtimes(L,T_shuffle,shuffles)
    true_LT=compute_loss(L,T)
    return calculate_pvalue(true_LT,loss_list_LT)

def test_3(L,A,B,shuffles):
    '''
    this function conducts the third test, returns the p-value and the overlap score 
    '''
    A_shuffle=np.copy(A)
    B_shuffle=np.copy(B)
    loss_list_Bresidual=stratify_B_n_times_diff(L,A_shuffle,B_shuffle,shuffles,algo)
    #calculate the FI score for the original sample 
    y_pred_ones,y_pred_zeros=compute_third_testloss(L,A,B,algo)
    diff=get_prob(L,A)
    true_LBresidual=(sum(np.multiply(abs(y_pred_zeros-y_pred_ones),diff))/len(y_pred_ones))
    cutoff=sum(diff)/len(diff)
    return [calculate_pvalue(true_LBresidual,loss_list_Bresidual),cutoff]


def main_call(i,child_seed):
    '''
    this function is called for each of the samples using Pooling, returns the sample index,p-values and overlap score
    '''
    rng=np.random.default_rng(child_seed)
    sample_seed=rng.integers(2**32 - 1)
    np.random.seed(sample_seed)
    random.seed(sample_seed)
    A=np.array(dataset_linear[i][1])
    B=np.array(dataset_linear[i][2])
    L=np.array(dataset_linear[i][0])
    LA_p=test_1(L,A,shuffles)
    LB_p=test_1(L,B,shuffles)
    AB_p,cutoff=test_3(L,A,B,shuffles)
    return [i,LA_p,LB_p,AB_p,cutoff]


if __name__ == '__main__':
    #generate a seed and save it to the start of the output file 
    if(len(sys.argv)==6): #it means seed is already provided 
        ss = int(sys.argv[5])
    else:
	ss = SeedSequence()
    f=open(outputname,"a")
    f.write("Seed "+ str(ss.entropy))
    f.write("\n")
    f.close()
    #unique child seeds for each of the sample 
    child_seeds = ss.spawn(len(dataset_linear))
    #parallelizing , here the number of workers is set as default to the number of cpus, you can modify it  
    with Pool( initargs=(dataset_linear,)) as pool:
        res=pool.starmap(main_call, zip(range(len(dataset_linear)), child_seeds))
    
    df=pd.DataFrame(res)
    
    et = time.time()
    elapsed_time = et - st
    print("Pooling of "+outputname)
    print('start time: ', st,' end time: ',et, 'Execution time:', elapsed_time, 'seconds')
    df.to_csv(outputname, header=False, index=False,mode='a')
    


