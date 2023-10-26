import os 
os.environ['OMP_NUM_THREADS'] = '1'  ## This will prevent KRR algo from thrashing the server
from scipy.stats import norm
import random
import numpy as np
from sklearn.svm import SVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.metrics import mean_squared_error
from sklearn.neural_network import MLPRegressor
from sklearn.utils._testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
@ignore_warnings(category=ConvergenceWarning)
def normalize(x_test,y_test):
    '''
    Function to standard normalize A|L for all the values of L to mean of 0 and variance 1 
    '''
    unique_values=np.unique(x_test)
    for value in unique_values:
        indices = np.where(x_test == value)[0]
        sigma = np.std(y_test[indices])
        y_test[indices] = (y_test[indices])/sigma
    return y_test 
    

def nlr_train_predict(xG, yG, algo, xL=None):
    '''
    Non-linear regression (NLR); 
    uses input gene xG and optionally input genotype xL to predict output gene yG, using the specified NLR algo
    returns predicted values and trained model (regressor)
    '''
    if xL is None :
        x = xG.reshape(-1,1)
    else:
        x = np.column_stack((xL, xG))
        
            
    if(algo=="SVR"):  #SVM
        regressor = SVR(kernel = 'rbf')
    elif(algo=="KRR"):
        regressor = KernelRidge(kernel = 'rbf')
    elif(algo=="ANN"):
        regressor=MLPRegressor()
    else:
        print("Invalid Algorithm")
        assert False
        
    regressor.fit(x,yG)
    y_pred = regressor.predict(x)
        
    return y_pred,regressor


def calculate_pvalue(original, loss_list):
    '''
    calculate the p value (with +1/+2 correction) 
    '''
    pvalNr = sum(i <= original for i in loss_list)
    return (pvalNr + 1)/(len(loss_list) + 2)
    

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
    diff=[]
    for i in range(len(probs)):
        for j in range(i+1,len(probs)):
            diff.append(np.minimum(probs[i],probs[j])/norm_const)

    return diff


def FI_score(x,y,overlap ):
    '''
    Function to calcuate the overlap score 
    '''

    # Iterate over the pairs of lists and calculate element-wise absolute differences
    diff = np.abs(np.subtract(x, y))

    return np.sum(np.multiply(diff,overlap))/len(x)


def compute_Luniqs_predns(L,A,B,algo):
    '''
    Function to calculate the predictions of y for different values of L in the 4th (or 2nd) test 
    '''
    y_pred = []
    _, regressor = nlr_train_predict(A, B, algo, L)
    
    
    unique_values = np.unique(L)
    assert (len(unique_values)==2 and (unique_values==[0,1]).all()) or (len(unique_values)==3 and (unique_values==[0,1,2]).all())
    for value in unique_values:
        L_1 = np.full_like(L, value)
        X_ = np.column_stack((L_1, A))
        y_predict_ = regressor.predict(X_)
        y_pred.append(y_predict_)
    
    return y_pred


def stratify_permute_variable(L, variable):
    '''
    Function to stratify permute a variable 
    '''
    unique_values = np.unique(L)
    permuted_variable = np.empty_like(variable)

    for value in unique_values:
        indices = np.where(L == value)[0]
        permuted_indices = np.random.choice(indices, len(indices),replace=False)
        permuted_variable[indices] = variable[permuted_indices]

    return permuted_variable


def test_4(L,A,B,shuffles,algo):
    '''
    Function for the 4th test , same function is used for the second test 
    '''

    overlap=get_prob(L,A)
    perm_loss=[]
    for i in range(shuffles):
            y_pred=compute_Luniqs_predns(L,A,stratify_permute_variable(L,B),algo)
            
            total_FI=0
            count=0
            for i in range(len(y_pred)):
                for j in range(i + 1, len(y_pred)):
                    total_FI+=FI_score(y_pred[i],y_pred[j],overlap[count])
                    count+=1
            assert count>0
            total_FI/=count
            perm_loss.append(total_FI)
    
    y_pred_original=compute_Luniqs_predns(L,A,B,algo)
    assert len(y_pred)==len(y_pred_original)
    assert len(y_pred)==2 or len(y_pred)==3
    original_loss=0
    count=0
    for i in range(len(y_pred_original)):
        for j in range(i + 1, len(y_pred_original)):
            original_loss+=FI_score(y_pred_original[i],y_pred_original[j],overlap[count])
            count+=1
    assert count>0
    original_loss/=count
    overlap_scores = [np.sum(overlaps) / len(overlaps) for overlaps in overlap]
    total_overlap_score = np.sum(overlap_scores)/count
    return [calculate_pvalue(original_loss,perm_loss),total_overlap_score,original_loss,perm_loss]


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
    #return NLL loss 
    #return np.mean(((((y_test-y_mean)/y_std)**2)/2)+np.log(2*np.pi)+np.log((y_std)))
    return -np.mean(np.log(np.exp(-(((y_test-y_mean)/y_std)**2)/2)/np.sqrt(2*np.pi*(y_std**2))))


def test_1(L,B,shuffles):
    '''
    Function for the first test 
    '''
    perm_loss=[]
    for i in range(0,shuffles):
        perm_loss.append(compute_1_loss(L,np.random.permutation(B)))

    original_loss=compute_1_loss(L,B)
    return calculate_pvalue(original_loss,perm_loss),original_loss,perm_loss    


def test_3_loss(A,B,algo):
    '''
    function to return the loss in the third function 
    '''
    Bpred, _ = nlr_train_predict(A, B, algo)
    mse = mean_squared_error(B, Bpred)
    return mse

    
def test_3(L,A,B,shuffles,algo):
    '''
    function for the third test
    '''
    unique_values = np.unique(L)
    p_values = []
    original_losses=[]
    perm_losslist=[]
    for value in unique_values:
        indices = np.where(L == value)[0]
        A_value = A[indices]
        B_value = B[indices]
        
        original_loss = test_3_loss(A_value, B_value, algo)
        original_losses.append(original_loss)
        perm_losses = [test_3_loss(A_value, np.random.permutation(B_value), algo) for _ in range(shuffles)]
        perm_losslist.append(perm_losses)
        p_value = calculate_pvalue(original_loss, perm_losses)
        p_values.append(p_value)
    
    #min_original_loss=original_losses[p_values.index(min(p_values))]    
    min_p_value = min(p_values)
    return min_p_value,original_losses,perm_losslist


def test_2(L,A,B,shuffles,algo):
    Apred, _ = nlr_train_predict(B, A, algo)
    Aresid = A - Apred
    out,original_loss,perm_loss = test_1(L, Aresid, shuffles)
    return out,original_loss,perm_loss    


def combine_tests(L,A,B,shuffles,algo,normal):
    '''
    Function to combine all the tests 
    '''
    if( len(np.unique(L))==1): # fix the number 10
        print("Zero variance of L, skipping trio")
        return [None]
    if( 1 in [sum(L==x) for x in np.unique(L)] ):
        print(" Only single value for a genotype value, can't do the statistics ")
        return [None]
    if(normal==True):
        A=normalize(L,A)
        B=normalize(L,B)

    T1_p,T1_OL,T1_PL=test_1(L,B,shuffles)
    T2_p,T2_OL,T2_PL=test_2(L,A,B,shuffles,algo)
    T3_p,T3_OL,T3_PL=test_3(L,A,B,shuffles,algo)
    T4_p,T4_OS,T4_OL,T4_PL=test_4(L,A,B,shuffles,algo)
    p_final=np.max([T1_p,T2_p,T3_p,T4_p])
    return [p_final,T1_p,T2_p,T3_p,T4_p,T4_OS,T1_OL,T2_OL,T3_OL,T4_OL,T1_PL,T2_PL,T3_PL,T4_PL]
