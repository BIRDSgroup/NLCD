# %%
import matplotlib.pyplot as plt
import seaborn as sns
import random
import numpy as np
import pandas as pd
from sklearn.metrics import (precision_recall_curve,PrecisionRecallDisplay)
from sklearn.metrics import precision_score,recall_score
from sklearn import metrics 
import seaborn as sns
import os 
import sys

# %%
def calculate_precision_recall(res,ind,permutations=100,iscit=False):
    precision=[]
    recall=[]
    if iscit==True:
        #initial row is the seed value 
        citpresults= pd.read_csv(res,skiprows=[1])
        #set the gorund truth as 1
        citpresults.loc[:,'p_res'] = 1 
        citpresultsInd= pd.read_csv(ind,skiprows=[1])
        #set the ground truth as 0 
        citpresultsInd.loc[:,'p_res'] = 0
        pCitMix=pd.concat([citpresults,citpresultsInd],ignore_index=True)
        pCitMix.columns= ['pred','p_TL', 'p_TG', 'p_GL','p_Lind','GroundTruth']
        uniq_val=np.unique(pCitMix['pred'])
       
        for i in uniq_val:
            pCitMix['results']=pCitMix.apply (lambda row: 1 if row['pred'] <=i   else 0, axis=1)
            precision.append(precision_score(pCitMix['GroundTruth'],pCitMix['results']))
            recall.append(recall_score(pCitMix['GroundTruth'],pCitMix['results']))
    else:
        presults=pd.read_csv(res,header=None,skiprows=1)
        presults = presults.iloc[: , :-1]
        presults.loc[:,0]=1
        presultsInd=pd.read_csv(ind,header=None,skiprows=1)
        presultsInd = presultsInd.iloc[: , :-1]
        presultsInd.loc[:,0]=0
        pMix=pd.concat([presults,presultsInd],ignore_index=True)
        pMix.columns= ['GroundTruth','p_L->A', 'p_L->B', 'conditional']
        pMix['pred']=pMix.apply (lambda row: max(row[1],row[2],row[3]) , axis=1) 
        #correcting for p value of 0 and p value of 1 
        pMix['pred']=pMix.apply (lambda row: (row['pred']*permutations+1)/(permutations+2) , axis=1)
        uniq_val=np.unique(pMix['pred'])
        for i in uniq_val:
            pMix['results']=pMix.apply (lambda row: 1 if row['pred'] <=i   else 0, axis=1)
            precision.append(precision_score(pMix['GroundTruth'],pMix['results']))
            recall.append(recall_score(pMix['GroundTruth'],pMix['results']))
            
    return precision,recall 
        

color_algo={'SVR':(0.8941176470588235, 0.10196078431372557, 0.10980392156862737),'ANN':(0.21568627450980393, 0.4941176470588236, 0.7215686274509804),'CIT': (0.3019607843137256, 0.6862745098039216, 0.29019607843137263),'KRR':(0.5960784313725492, 0.3058823529411765, 0.6392156862745098)}
# %%
### yeast 
precision_nlmr_SVR,recall_nlmr_SVR= calculate_precision_recall("./../results/yeast/nlmr/causal1kp100.txt","./../results/yeast/nlmr/indep1kp100.txt")
precision_CIT,recall_CIT= calculate_precision_recall("./../results/yeast/cit/yeastcausal1kcit.csv","./../results/yeast/cit/yeastindp1kcit.csv",iscit=True)
fig, ax = plt.subplots()
disp=PrecisionRecallDisplay(precision=precision_CIT, recall=recall_CIT)
disp.plot(ax=ax,name="CIT",color=color_algo['CIT'])
disp=PrecisionRecallDisplay(precision=precision_nlmr_SVR, recall=recall_nlmr_SVR)
disp.plot(ax=ax, name="SVR", color=color_algo['SVR'])
ax.legend( loc="best")
ax.set_title("P-R curve for yeast data No of samples= 1000")
plt.savefig("./prcurveyeast1k.png")



