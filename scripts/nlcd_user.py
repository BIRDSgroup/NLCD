from nlcd_main import *
import time 
import sys
from numpy.random import SeedSequence
from multiprocessing import Pool
import pandas as pd
def read_data(data):

    #read the input file 
    fo=open(data, "r")
    L=[]
    A=[]
    B=[]
# the simulated data is generated with a seed at the start of the file
# so we skip it, but not for the yeast data 
# if yeast is present then dont read the seed line
    if(data.find('yeast')==-1 and data.find('human')==-1): #if yeast is not present
        line=fo.readline() # read the line, the control will start from the next line 

    while(1):
    # the first line of the trio is the parameter configuration, skip it 
        line=fo.readline()
    # read the genotype information 
        line=fo.readline()
        if line== "" :
            break
        l = [j for j in line.split()]
        L.append(np.array([int(i) for i in l]))
        line=fo.readline()
        a = [j for j in line.split()]
        A.append(np.array([float(i) for i in a]))
        line=fo.readline()
        b = [j for j in line.split()]
        B.append(np.array([float(i) for i in b]))
#combine all the samples into a list 
    dataset = [i for i in zip(L,A,B)]
    fo.close()
    return dataset
def read_configuration(data):
    fo=open(data, "r")
    confignames=[]
    if(data.find('yeast')==-1 and data.find('human')==-1): #if yeast is not present
        line=fo.readline() # read the line, the control will start from the next line 
    while(1):
        line=fo.readline()
        if(line==""):
            break
        line=line.replace("\n","")
        line=line.replace("\"","")
        line=line.split(" ")  
        configs=[]  
        for i in line: 
            configs.append(i)
        confignames.append(configs)
        line=fo.readline()
        line=fo.readline()
        line=fo.readline()
    fo.close()
    return pd.DataFrame(confignames)

def nlcd_batch(shared_data, shuffles, algo, reverse=False,sample_seed=None,normal=False,usercall=False):
    st = time.time()
    if(sample_seed==None):
        ss = SeedSequence()
    else:
        ss = SeedSequence(sample_seed)
    print("Seed Set at ",ss.entropy)
    
    # unique child seeds for each of the trio [cf. https://numpy.org/doc/stable/reference/random/parallel.html#seedsequence-spawning ]
    ntrios = len(shared_data)
    child_seeds = ss.spawn(ntrios)
    child_seeds_ints = []
    for i in range(ntrios):
        rng = np.random.default_rng(child_seeds[i])
        child_seeds_ints.append(rng.integers(sys.maxsize))  
    verbose = False    
    #parallelizing, here the number of workers is set as default to the number of cpus, you can modify it  
    with Pool() as pool:
        res = pool.starmap(nlcd_single_for_batch, zip(shared_data, [shuffles]*ntrios, [algo]*ntrios, child_seeds_ints, [verbose]*ntrios, [reverse]*ntrios,[normal]*ntrios))
    #print(res)
    data=[[sublist[0], sublist[1], sublist[2], sublist[3],sublist[4],sublist[5],sublist[14]] if sublist!=[None] else [None] for sublist in res]
    p_df=pd.DataFrame(data,dtype=object)
    data=[[sublist[6],*sublist[10]] if sublist!=[None] else [None] for sublist in res] # Test 1
    t1loss=pd.DataFrame(data).T
    data=[[sublist[7],*sublist[11]] if sublist!=[None] else [None] for sublist in res] # Test 2
    t2loss=pd.DataFrame(data).T
    data=[[sublist[9],*sublist[13]] if sublist!=[None] else [None] for sublist in res] # Test 4
    t4loss=pd.DataFrame(data).T
    data=[list([x, *y] for x, y in zip(sublist[8], sublist[12]))[0] if sublist!=[None] else [None] for sublist in res]
    t3loss_0=pd.DataFrame(data).T
    data=[list([x, *y] for x, y in zip(sublist[8], sublist[12]))[1] if sublist!=[None] else [None] for sublist in res]
    t3loss_1=pd.DataFrame(data).T
    data=[list([x, *y] for x, y in zip(sublist[8], sublist[12]))[2] if sublist!=[None] and len(list([x, *y] for x, y in zip(sublist[8], sublist[12])))>2   else [None] for sublist in res]
    t3loss_2=pd.DataFrame(data).T
    #df=pd.DataFrame(res)
    et = time.time()
    elapsed_time = et - st
    print("Algo ",algo," shuffles ",shuffles," datasize = ",len(shared_data[0][0])," reverse ",reverse)
    print('start time: ', st,' end time: ',et, 'Execution time:', elapsed_time, 'seconds')
    p_df.columns=['p_final','p_LassocB','p_LassocA|B','p_AassocB|L','p_LindB|A','OS Test 4','child_seed']
    p_df['parent_seed']=[ss.entropy]+['same']*(ntrios-1)
    if usercall==False:
        return p_df,t1loss,t2loss,t4loss,t3loss_0,t3loss_1,t3loss_2
    else:
        return p_df

def nlcd_single_for_batch(singletriodata, shuffles, algo, sample_seed=None, verbose=False, reverse=False, normal=False):
    assert len(singletriodata)==3
    return nlcd_single(singletriodata[0],singletriodata[1],singletriodata[2], shuffles, algo, sample_seed, verbose, reverse,normal)

def nlcd_single(L, A, B, shuffles, algo, sample_seed=None, verbose=True, reverse=False,normal=False):
    if(sample_seed==None):
        rng = np.random.default_rng()
        sample_seed = rng.integers(sys.maxsize) 
    rng=np.random.default_rng(sample_seed)
    
    sample_seed1 = rng.integers(2**32 - 1) # numpy will throw an error if the seed is set as sys.maxsize  
    
    np.random.seed(sample_seed1)
    ## check if L is haploid but contains only two unique values 
    if (2 in L and len(np.unique(L))==2) :
        print("Diploid but only 2 unique values for L hence treating it as haploid")
        L=np.array([0 if (x == 0 or x == 1) else 1 for x in L])
    if reverse==False:
        out=combine_tests(L,A,B,shuffles,algo,normal)
    elif reverse==True:
        out=combine_tests(L,B,A,shuffles,algo,normal)
    else:
        print("Invalid entry for reverse parameter")
    if out!=[None]:
        out.append(sample_seed)
    if verbose==True and out!=[None]:
        print("The final p value is ",out[0])
        print("Test 1 L assoc B ",out[1])
        print("Test 2 L assoc A | B ",out[2])
        print("Test 3 A assoc B | L ",out[3])
        print("Test 4 L ind B | A ",out[4])
        print("Overlap score from Test 4 ",out[5])
        print("Seed set at ",out[14])
        return out[0:6]+[out[14]]
    elif verbose==True:
        print("Output vector is None")
        return 
    
    return out


