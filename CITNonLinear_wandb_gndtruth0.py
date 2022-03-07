#!/usr/bin/env python
# coding: utf-8

# In[19]:


yeast=open("yeast_residual_data_full_1000_gt_2.txt","r")



#yeast data read 
L=[]
A=[]
B=[]
for i in range(0,1000):
  line=yeast.readline()
  #line=line[1:-2] #remove double quotes 
  #param = [j for j in line.split()]
  #print(param)
  #chrname.append(param[1])
  #g1.append(param[2])
  #g2.append(param[3])
  line=yeast.readline()
  l = [j for j in line.split()]
  L.append([int(i) for i in l])
  line=yeast.readline()
  a = [j for j in line.split()]
  A.append([float(i) for i in a])
  line=yeast.readline()
  b = [j for j in line.split()]
  B.append([float(i) for i in b])
dataset_0 = [i for i in zip(L,A,B)]


# In[20]:


#yeast data
yeast=open("yeast_residual_data_full_1000_gt_1.txt","r")



#yeast data read 
L=[]
A=[]
B=[]
for i in range(0,1000):
  line=yeast.readline()
  #line=line[1:-2] #remove double quotes 
  #param = [j for j in line.split()]
  #print(param)
  #chrname.append(param[1])
  #g1.append(param[2])
  #g2.append(param[3])
  line=yeast.readline()
  l = [j for j in line.split()]
  L.append([int(i) for i in l])
  line=yeast.readline()
  a = [j for j in line.split()]
  A.append([float(i) for i in a])
  line=yeast.readline()
  b = [j for j in line.split()]
  B.append([float(i) for i in b])
dataset_1 = [i for i in zip(L,A,B)]

#reshapefunction


# In[7]:



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
import wandb
from sklearn.feature_selection import mutual_info_regression
from scipy import stats


# In[8]:


from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
nlcor=importr('nlcor')
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()


# In[9]:


class MDN_module(tf.keras.Model):

    def __init__(self, neurons=15, components = 1):
        super(MDN_module, self).__init__(name="MDN_module")
        self.neurons = neurons
        self.components = components

        #chaging activation to relu from linear, changin relu to sigmoid 
        for i in range(1,3):
          s="self"+".h"+str(i)+"= Dense(neurons, activation=\"relu\", name="+"'h"+str(i)+"')"
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


# In[10]:


no_parameters=3
components=1

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
    gm = tfd.MixtureSameFamily(
        mixture_distribution=tfd.Categorical(probs=alpha),
        components_distribution=tfd.Normal(
            loc=mu,       
            scale=sigma))
    log_likelihood = gm.log_prob(tf.transpose(y))
    return -tf.reduce_mean(log_likelihood, axis=-1)


def eval_mdn_model(x_test, y_test, mdn_model):
    y_pred = mdn_model.predict(np.reshape(x_test,newshape=(len(x_test),-1)))
    alpha,mu,sigma = slice_parameter_vectors(y_pred)
    mdn_nll = gnll_eval(y_test.astype(np.float32),alpha, mu, sigma).numpy()
    return mdn_nll
#reshapefunction
def reshapevar(X):
  return np.reshape(X,newshape=(len(X),-1))


# In[11]:


def compute_loss(P,Q):
  opt = tf.optimizers.Adam(1e-2)
  mdn_PQ = MDN_module()
  mdn_PQ.compile(loss=gnll_loss, optimizer=opt)
  mdn_PQ.fit(x=reshapevar(P), y=np.array(Q).T,epochs=100,  batch_size=64,verbose=0)
  return eval_mdn_model(P,Q,mdn_PQ)


# In[12]:


def shuffleBtimes(P,Q,B):
  loss=[]
  for i in range(0,B):
    loss.append(compute_loss(P,np.random.permutation(Q)))
  return loss


# In[13]:


def LinearLABData():
  L = np.random.binomial(1,0.5,1000)  
  beta0 = np.ones(1000)-0.4
  #beta1 = 0.5
  beta1=3
  beta2= 0.3
  beta3=0.8
  eps0 = np.random.standard_normal(1000)
  eps1 = np.random.standard_normal(1000)
  A = beta0 + beta1*L + eps0
  B = beta2+ beta3*A+ eps1 
  #plt.scatter(A,B)
  #plt.title("A vs B")
  #plt.xlabel("A")
  #plt.ylabel("B")
  return [L,A,B]


# In[14]:


#not using it now 
def residual(P,Q):
  model=sm.OLS(Q,P).fit()
  return model.resid
  


# In[15]:


#yeast_name=""
def yeast_data(i,ind):
    global yeast_name
    yeast_name="yeast_"+str(i)+"_"+str(ind)
    ds = eval("dataset_"+str(i)+"["+str(ind)+"]")
    L_dist = np.array(ds[0]) #np.array(ds[0])
    A_dist = np.array(ds[1])
    B_dist = np.array(ds[2])
    #plt.scatter(A_dist,B_dist)
    #plt.title("A vs B")
    #plt.xlabel("A")
    #plt.ylabel("B")
    return [L_dist,A_dist,B_dist]


# In[16]:


def calculate_pvalue(original,loss_list):
    return sum(i < original for i in loss_list)/len(loss_list)
    


# In[17]:


wandb.init(
      name="yeast_simplenn_groundtruth_0_100shuffles",
      project="yeast_pvalue_otherstats",
      entity="birds")


# In[24]:


random.seed(24)
shuffles=100
for i in range(0,1000):
#L,A,B=LinearLABData()
    L,A,B=yeast_data(0,i)
    B_resid=residual(A,B)
    A_shuffle=np.copy(A)
    B_shuffle=np.copy(B)
    B_resid_shuffle=np.copy(B_resid)
    loss_list_LA=shuffleBtimes(L,A_shuffle,shuffles)
    loss_list_LB=shuffleBtimes(L,B_shuffle,shuffles)
    #loss_list_LindB_A=shuffleBtimes(L,B_resid_shuffle,shuffles) #conditional independence test
    true_LA=compute_loss(L,A)
    true_LB=compute_loss(L,B)
    #print(true_LB)
    #print(loss_list_LB)
    loss_LA_p=calculate_pvalue(true_LA,loss_list_LA)
    loss_LB_p=calculate_pvalue(true_LB,loss_list_LB)
    loss_LA_mutual=mutual_info_regression(L.reshape(-1,1),A)
    loss_LB_mutual=mutual_info_regression(L.reshape(-1,1),B)
    loss_LA_spear=stats.spearmanr(L,A)[0]
    loss_LB_spear=stats.spearmanr(L,B)[0]
    loss_LA_pear=stats.pearsonr(L,A)[0]
    loss_LB_pear=stats.pearsonr(L,B)[0]
    loss_LA_nlcor=np.array(nlcor.nlcor(L,A)[0])[0]
    loss_LB_nlcor=np.array(nlcor.nlcor(L,B)[0])[0]
    
    #true_LindB_A=compute_loss(L,B_resid)
    wandb.log({"dataset":i,"loss_LA_p":loss_LA_p,"loss_LB_p":loss_LB_p,"loss_LA_mutual":loss_LA_mutual,"loss_LB_mutual":loss_LB_mutual,"loss_LA_spear":loss_LA_spear,"loss_LB_spear":loss_LB_spear,"loss_LA_pear":loss_LA_pear,"loss_LB_pear":loss_LB_pear,"loss_LA_nlcor":loss_LA_nlcor,"loss_LB_nlcor":loss_LB_nlcor,"true_LA":true_LA,"true_LB":true_LB})


# In[25]:


wandb.finish()

