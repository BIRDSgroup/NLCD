#!/usr/bin/env python
# coding: utf-8

# In[9]:


#loading yeast groundtruth 0 data
yeast=open("../yeast_residual_data_full_209k_gt2.txt","r")



#yeast data read 
L=[]
A=[]
B=[]
#209157,62296
for i in range(0,209157):
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


# In[18]:



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
import pandas as pd 


# In[5]:


from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
nlcor=importr('nlcor')
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()


# In[44]:


df = pd.DataFrame(columns = ['Dataset', 'pearson','spearman','Mututal_info','nlcor'])


# In[48]:


for i in range(len(dataset_0)):
    A = np.array(dataset_0[i][1])
    B = np.array(dataset_0[i][2])
    loss_mutual=mutual_info_regression(A.reshape(-1,1),B)[0]
    loss_spear=stats.spearmanr(A,B)[0]
    loss_pear=stats.pearsonr(A,B)[0]
    #print(loss_mutual,loss_spear,loss_pear)
    loss_nlcor=np.array(nlcor.nlcor(A,B)[0])[0]
    df.loc[i] = [i]+ [loss_pear]+[loss_spear]+[loss_mutual]+[loss_nlcor]
    


# In[52]:


df.to_csv("Gt2_full_corr.csv")

