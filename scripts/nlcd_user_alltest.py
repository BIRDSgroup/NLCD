from nlcd_main import *
import time
import sys
from numpy.random import SeedSequence
from multiprocessing import Pool
import pandas as pd 
def combine_tests(L,A,B,shuffles,algo):
    '''
    Function to combine all the tests 
    '''
    #LA_p=test_1(L,A,shuffles)
    LB_p=test_1(L,B,shuffles)
    LAgvnBv1=test_2(L,A,B,shuffles,algo,version=1)
    LAgvnBv2,overlapscoret2v2=test_2(L,B,A,shuffles,algo,version=2) # v2
    LAgvnBv3,overlapscoret2v3=test_2(L,A,B,shuffles,algo,version=3)
    ABgvnL=test_3(L,A,B,shuffles,algo)
    LindBgvnA,overlapscoret4=test_4(L,A,B,shuffles,algo)
    return [LB_p,LAgvnBv1,LAgvnBv2,LAgvnBv3,ABgvnL,LindBgvnA,overlapscoret2v2,overlapscoret2v3,overlapscoret4]

def nlcd_batch(shared_data, shuffles, algo, reverse=False,sample_seed=None):
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
        res = pool.starmap(nlcd_single_for_batch, zip(shared_data, [shuffles]*ntrios, [algo]*ntrios, child_seeds_ints, [verbose]*ntrios, [reverse]*ntrios))
    
    df=pd.DataFrame(res)
    et = time.time()
    elapsed_time = et - st
    print('start time: ', st,' end time: ',et, 'Execution time:', elapsed_time, 'seconds')
    df.columns=['LB_p','LA|Bv1','LA|Bv2','LA|Bv3','AB|L','LindB|A','OSt2v2','OSt2v3','OSt4','seed']
    #save the seed in the end 
    df['parent_seed']=[ss.entropy]+['same']*(ntrios-1)
    return df

def nlcd_single_for_batch(singletriodata, shuffles, algo, sample_seed=None, verbose=False, reverse=False):
    assert len(singletriodata)==3
    return nlcd_single(singletriodata[0],singletriodata[1],singletriodata[2], shuffles, algo, sample_seed, verbose, reverse)

def nlcd_single(L, A, B, shuffles, algo, sample_seed=None, verbose=True, reverse=False):
    if(sample_seed==None):
        rng = np.random.default_rng()
        sample_seed = rng.integers(sys.maxsize) 
    rng=np.random.default_rng(sample_seed)
    
    sample_seed1 = rng.integers(2**32 - 1) # numpy will throw an error if the seed is set as sys.maxsize  
    np.random.seed(sample_seed1)
    ## check if L is haploid but contains only two unique values 
    if (2 in L and len(np.unique(L))==2) :
        print("Hapoloid but only 2 unique values for L hence treating it as diploid")
        L=np.array([0 if (x == 0 or x == 1) else 1 for x in L])
    if reverse==False:
        out=combine_tests(L,A,B,shuffles,algo)
    elif reverse==True:
        out=combine_tests(L,B,A,shuffles,algo)
    else:
        print("Invalid entry for reverse parameter")
    out.append(sample_seed)
    if verbose==True:
	#need to modify these print statements 
        print("The final p value is ",out[0])
        print("Test 1 L assoc B ",out[1])
        print("Test 2 L assoc A | B ",out[2])
        print("Test 3 A assoc B | L ",out[3])
        print("Test 4 L ind B | A ",out[4])
        print("Overlap score from Test 2 ",out[5])
        print("Overlap score from Test 4 ",out[6])
        print("Seed set at ",out[7])
        
    
    return out


