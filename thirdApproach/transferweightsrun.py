#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#loading the libraries
from tensorflow_probability import distributions as tfd
import tensorflow as tf
from tensorflow.keras.layers import Input, Dense, Activation, Concatenate

from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from sklearn.datasets import load_boston
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns
import random
import numpy as np
import statsmodels.api as sm
import pickle
from sklearn.feature_selection import mutual_info_regression
from scipy import stats
import rpy2


# In[ ]:


f=open("Linearbest4.txt","a")


# In[ ]:


import sys
#First argument for the number of runs, second for the initial run
j=int(sys.argv[2])
n=int(sys.argv[1])


# In[ ]:


#defining the class MDN
class MDN_module(tf.keras.Model):

    def __init__(self, neurons=15, components = 1):
        super(MDN_module, self).__init__(name="MDN_module")
        self.neurons = neurons
        self.components = components

        #chaging activation to relu from linear, changin relu to sigmoid 
        for i in range(1,3):
          s="self"+".h"+str(i)+"= Dense(neurons, activation=\"relu\",name="+"'h"+str(i)+"')"
          exec(s)
        self.alphas = Dense(components, activation="softmax", name="alphas")
        self.mus = Dense(components, activation="linear",name="mus") 
        self.sigmas = Dense(components, activation="nnelu",name="sigmas") #activation changed from linear to default
        self.pvec = Concatenate(name="pvec")
        
    def call(self, inputs):
        x=self.h1(inputs)
        #x=self.inputA(inputs)
        x=self.h2(x)
        alpha_v = self.alphas(x)
        mu_v = self.mus(x)
        sigma_v = self.sigmas(x)
        
        return self.pvec([alpha_v,mu_v, sigma_v])


# In[ ]:


no_parameters=3
components=1

epochs=100
batchsize=32
def nnelu(input):
    """ Computes the Non-Negative Exponential Linear Unit
    """
    return tf.add(tf.constant(1, dtype=tf.float32), tf.nn.elu(input))

def slice_parameter_vectors(parameter_vector):
    """ Returns an unpacked list of paramter vectors.
    """
    return [parameter_vector[:,i*components:(i+1
    )*components] for i in range(no_parameters)]

def gnll_loss(y, parameter_vector):
    """ Computes the mean negative log-likelihood loss of y given the mixture parameters.
    """
    alpha,mu,sigma = slice_parameter_vectors(parameter_vector) # Unpack parameter vectors
    #tf.print(sigma)
    gm = tfd.MixtureSameFamily(
           mixture_distribution=tfd.Categorical(probs=alpha),
           components_distribution=tfd.Normal(
           loc=mu,       
           scale=sigma))
    
    
    
    log_likelihood =  gm.log_prob(tf.transpose(y)) # Evaluate log-probability of y 
    return -tf.reduce_mean(log_likelihood, axis=-1) 

tf.keras.utils.get_custom_objects().update({'nnelu': Activation(nnelu)})

def gnll_eval(y,alpha, mu, sigma):
    """ Computes the mean negative log-likelihood loss of y given the mixture parameters.
    """
    #print(alpha)
    gm = tfd.MixtureSameFamily(
        mixture_distribution=tfd.Categorical(probs=alpha),
        components_distribution=tfd.Normal(
            loc=mu,       
            scale=sigma))
    log_likelihood = gm.log_prob(tf.transpose(y))
    return -tf.reduce_mean(log_likelihood, axis=-1)


def eval_mdn_model(x_test, y_test, mdn_model):
    """ Evaluate the model to get the loss for the given x and y 
    """
    y_pred = mdn_model.predict(np.reshape(x_test,newshape=(len(x_test),-1)))
    alpha,mu,sigma = slice_parameter_vectors(y_pred)
    mdn_nll = gnll_eval(y_test.astype(np.float32),alpha, mu, sigma).numpy()
    return mdn_nll
#reshapefunction
def eval_mdn_model_mle(x_test,y_test):
        indices_1 = [i for i, x in enumerate(x_test) if x == 1]
        indices_0 = [i for i, x in enumerate(x_test) if x == 0]
        mu_0=np.mean(y_test[indices_0])
        mu_1=np.mean(y_test[indices_1])
        sigma_0=np.std(y_test[indices_0])
        sigma_1=np.std(y_test[indices_1])
        y_mean=np.zeros((len(y_test),1))
        y_mean[indices_1]=mu_1
        y_mean[indices_0]=mu_0
        y_std=np.zeros((len(y_test),1))
        y_std[indices_1]=sigma_1
        y_std[indices_0]=sigma_0
        alpha=np.ones((len(y_mean),1))
        return gnll_eval(y_test,alpha,y_mean,y_std).numpy()
    
def reshapevar(X):
  """
  Function to reshape the vector for the input 
  """
  return np.reshape(X,newshape=(len(X),-1))


# In[ ]:


def compute_loss(P,Q,mle=False):
    """ Compute the loss for the given pair
    """
    if(mle==False):
        opt = tf.optimizers.Adam(1e-2)
        mdn_PQ = MDN_module()
        mdn_PQ.compile(loss=gnll_loss, optimizer=opt)
        mdn_PQ.fit(x=reshapevar(P), y=np.array(Q).T,epochs=100,  batch_size=64,verbose=0)
        #return np.array(nlcor.nlcor(P,Q)[0])[0]
        return eval_mdn_model(P,Q,mdn_PQ)
    else:
        return eval_mdn_model_mle(P,Q)
        
        
        


# In[ ]:


def compute_loss_y_pred(P,Q,mle=False):
    """ Compute the loss for the given pair
    """
    if(mle==False):
        opt = tf.optimizers.Adam(1e-2)
        mdn_PQ = MDN_module()
        mdn_PQ.compile(loss=gnll_loss, optimizer=opt)
        mdn_PQ.fit(x=reshapevar(P), y=np.array(Q).T,epochs=100,  batch_size=64,verbose=0)
        y_pred = mdn_PQ.predict(np.reshape(P,newshape=(len(P),-1)))
        return y_pred[:,1]
    else:
        indices_1 = [i for i, x in enumerate(P) if x == 1]
        indices_0 = [i for i, x in enumerate(P) if x == 0]
        mu_0=np.mean(Q[indices_0])
        mu_1=np.mean(Q[indices_1])
        #sigma_0=np.std(Q[indices_0])
        #sigma_1=np.std(Q[indices_1])
        y_mean=np.zeros((len(Q),1))
        y_mean[indices_1]=mu_1
        y_mean[indices_0]=mu_0
        return y_mean.reshape((len(y_mean),))
        


# In[ ]:


def shuffleBtimes(P,Q,B,mle=False):
    """ Shuffle Q B times and compute the loss 
    """
    loss=[]
    if(mle==False):
        for i in range(0,B):
          loss.append(compute_loss(P,np.random.permutation(Q)))
    else:
        for i in range(0,B):
          loss.append(compute_loss(P,np.random.permutation(Q),True))
    return loss


# In[ ]:


def calculate_pvalue(original,loss_list):
    '''
    calculate the p value 
    '''
    return sum(i < original for i in loss_list)/len(loss_list)
    
    


# In[ ]:


fo=open("/data/users/cs20s037/CITNonLinear/LinearDifferentvalues/testing_writingvalues_Linear0to1.txt", "r")
L=[]
A=[]
B=[]
#fe=open("dataset_params.txt",'w')
for i in range(0,121):
    line=fo.readline()
    #fe.write(line)
    #line=line[1:-2] #remove double quotes 
    #param = [j for j in line.split()]
    #print(param)
    #chrname.append(param[1])
    #g1.append(param[2])
    #g2.append(param[3])
    line=fo.readline()
    l = [j for j in line.split()]
    L.append([int(i) for i in l])
    line=fo.readline()
    a = [j for j in line.split()]
    A.append([float(i) for i in a])
    line=fo.readline()
    b = [j for j in line.split()]
    B.append([float(i) for i in b])
dataset_linear = [i for i in zip(L,A,B)]
fo.close()
#fe.close()


# In[ ]:


#function to train the model given the inputs
def modeltrain(x,y,weights=None):
    opt = tf.optimizers.SGD(1e-3)
    
    mdn_PQ = MDN_module()
    mdn_PQ.compile(loss=gnll_loss, optimizer=opt)
    if(weights!=None):
        mdn_PQ.fit(x, y,epochs=1,batch_size = batchsize,verbose=0)
        for i in range(0,6):
            mdn_PQ.layers[i].set_weights(weights[i])
            
    #calculate the prediction and the loss 
    mdn_PQ.fit(x, y,epochs=epochs,batch_size=batchsize ,verbose=0)
    y_pred = mdn_PQ.predict(x)
    alpha,mu,sigma = slice_parameter_vectors(y_pred)
    gm = tfd.MixtureSameFamily(
            mixture_distribution=tfd.Categorical(probs=alpha),
            components_distribution=tfd.Normal(
                loc=mu,       
                scale=sigma))
    log_likelihood = gm.log_prob(y).numpy()
    loss_o = -tf.reduce_mean(log_likelihood, axis=-1).numpy()
    return mdn_PQ,loss_o


# In[ ]:


#function to calculate the weights 
def calculateweights(mdn_PQ):
    weights=[]
    weightofa=mdn_PQ.layers[0].get_weights()[0][0]
    #get the bias
    biasweight=mdn_PQ.layers[0].get_weights()[1]
    #set the weights of L as zero 
    weightofL=np.zeros(15,dtype=np.float32)
    #combine the weights of L and A
    combinedweight=np.array([weightofL,weightofa])
    #combine the weights and biases
    layer0weightrans=[combinedweight,biasweight]
    weights.append(layer0weightrans)
    for i in range(1,6):
        a=mdn_PQ.layers[i].get_weights()
        weights.append(a)
    return weights


# In[ ]:


#function to stratify the data
def stratifydata(L,B):
    indices_1 = [i for i, x in enumerate(L) if x == 1]
    #changin x==0 to x=-1
    indices_0 = [i for i, x in enumerate(L) if x == 0]
    B_dist_temp=np.zeros(len(B))
    mod_indices_1=random.sample(indices_1,len(indices_1))
    for i in range(len(indices_1)):
        B_dist_temp[indices_1[i]]=B[mod_indices_1[i]]

    mod_indices_0=random.sample(indices_0,len(indices_0))
    for i in range(len(indices_0)):
        B_dist_temp[indices_0[i]]=B[mod_indices_0[i]]
    return B_dist_temp
    


# In[ ]:


#function to split the data
def splitdata(L,A,B):
    Aalone_L1=[]
    Aalone_L0=[]
    Balone_L1=[]
    Balone_L0=[]
    for i in range(len(L)):
        if(L[i]==1):
            Aalone_L1.append(A[i])
            Balone_L1.append(B[i])
        else:
            Aalone_L0.append(A[i])
            Balone_L0.append(B[i])
    
    return [np.array(Aalone_L1),np.array(Aalone_L0),np.array(Balone_L0),np.array(Balone_L1)]


# In[ ]:


#checking one cycle 
def thirdtestonecycle(L,A,B):
    loss_Aalone=[]
    loss_AL=[]
    halfdata= splitdata(L,A,B)
    trainthirdtest(L,A,B,halfdata,loss_Aalone,loss_AL)
    trainthirdtest(L,A,B,halfdata,loss_Aalone,loss_AL,True,half=0)
    trainthirdtest(L,A,B,halfdata,loss_Aalone,loss_AL,True,half=1)
    trainthirdtest(L,A,B,halfdata,loss_Aalone,loss_AL,True)
    #print("Loss with A alone: Loss A",loss_Aalone)
    #print("Loss with A and L: Loss LA",loss_AL)
    a= 1 if min(loss_Aalone)-min(loss_AL) >=0 else 0
    return min(loss_Aalone)-min(loss_AL),a


# In[ ]:


def thirdtestoriginal(L,A,B):
    #print("Original data")
    #print("Loss Order: Random,L=0 tranfer, L= 1 tranfer, Full data transfer")
    loss,neg=thirdtestonecycle(L,A,B)
    #print(" Difference(minimum from Loss A - minimum from Loss LA) ",loss)
    return loss


# In[ ]:


def thirdtestperm(L,A,B,n):
    loss_perm=[]
    for i in range(0,n):
        #print("Permutation ",i)
        #print("Loss Order: Random,L=0 tranfer, L= 1 tranfer, Full data transfer")
        B_perm=stratifydata(L,B)
        loss,neg=thirdtestonecycle(L,A,B_perm)
        #print(" Difference(minimum from Loss A - minimum from Loss LA) ",loss)
        loss_perm.append(loss)
    return loss_perm


# In[ ]:


#function to train the data for each configuration 
def trainthirdtest(L,A,B,halfdata,loss_Aalone,loss_AL,transfer=False,half=None):
    AL=np.concatenate([L.reshape(-1,1),A.reshape(-1,1)],axis=1)
    if half==None:
        mdn,loss=modeltrain(reshapevar(A),B.T)
        loss_Aalone.append(loss)
        weights=None
        if transfer==True :
            weights=calculateweights(mdn)
        mdn,loss=modeltrain(AL,B.T,weights)
        loss_AL.append(loss)
    elif half==1:
        Aalone_L1=halfdata[0]
        Balone_L1=halfdata[3]
        mdn,loss=modeltrain(reshapevar(Aalone_L1),Balone_L1.T)
        loss_Aalone.append(loss)
        weights=None
        if transfer==True :
            weights=calculateweights(mdn)
        mdn,loss=modeltrain(AL,B.T,weights)
        loss_AL.append(loss)
    elif half==0:
        Aalone_L0=halfdata[1]
        Balone_L0=halfdata[2]
        mdn,loss=modeltrain(reshapevar(Aalone_L0),Balone_L0.T)
        loss_Aalone.append(loss)
        weights=None
        if transfer==True :
            weights=calculateweights(mdn)
        mdn,loss=modeltrain(AL,B.T,weights)
        loss_AL.append(loss)
    


# In[ ]:


for i in range(j,j+n): 
    A=np.array(dataset_linear[i][1])
    B=np.array(dataset_linear[i][2])
    L=np.array(dataset_linear[i][0])
    shuffles=100
    A_shuffle=np.copy(A)
    B_shuffle=np.copy(B)
    #print("Original",B_shuffle)
    loss_list_LA=shuffleBtimes(L,A_shuffle,shuffles,True)
    loss_list_LB=shuffleBtimes(L,B_shuffle,shuffles,True)
    #loss_list_Bresidual=stratify_B_n_times_diff(L,A_shuffle,B_shuffle,shuffles) #conditional independence test
    true_LA=compute_loss(L,A,True)
    true_LB=compute_loss(L,B,True)
    #true_LBresidual=calculate_difference(L,A,B)
    #loss_list_Bresidual,true_LBresidual=calculateLshuffle(L,A,B,shuffles)
    true_LBresidual=thirdtestoriginal(L,A,B)
    loss_list_Bresidual=thirdtestperm(L,A,B,shuffles)
    LA_p=calculate_pvalue(true_LA,loss_list_LA)
    LB_p=calculate_pvalue(true_LB,loss_list_LB)
    AB_p=calculate_pvalue(true_LBresidual,loss_list_Bresidual)
    f.write(str(i)+","+str(LA_p)+","+str(LB_p)+","+str(AB_p)+","+"\n")
    #pickle_items=[loss_list_LA,loss_list_LB,loss_list_Bresidual,true_LA,true_LB,true_LBresidual,LA_p,LB_p,AB_p]
    #file_name=str(dataset_names[i])+".pkl"
    #open_file = open("./DLresultspickle1000shuffle/"+file_name, "wb")
    #pickle.dump(pickle_items, open_file)
    #open_file.close()


# In[ ]:


f.close()
