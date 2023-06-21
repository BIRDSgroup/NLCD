from nlcd_main import *
import time 
def init_worker(shared_data):
    # declare scope of a new global variable
    global data
    # store argument in the global variable for this process
    data = shared_data
 
def main_call(i,child_seed,data,shuffles,algo):
    '''
    this function is called for each of the samples using Pooling, returns the sample index,p-values and overlap score
    '''
    rng=np.random.default_rng(child_seed)
    sample_seed=rng.integers(2**32 - 1)
    np.random.seed(sample_seed)
    random.seed(sample_seed)
    A=np.array(data[i][1])
    B=np.array(data[i][2])
    L=np.array(data[i][0])
    return combine_tests(L,A,B,shuffles,algo)

def nlcd_batch(shared_data, shuffles, algo, sample_seed=None):
    st = time.time()
    if(sample_seed==None):
        ss = SeedSequence()
    else:
        ss = SeedSequence(sample_seed)
    print("Seed Set at ",ss.entropy)
    
    # unique child seeds for each of the trio [cf. https://numpy.org/doc/stable/reference/random/parallel.html#seedsequence-spawning ]
    ntrios = len(shared_data)
    child_seeds = ss.spawn(ntrios)
    child_seeds_ints = []; 
    for i in range(ntrios):
        rng = np.random.default_rng(child_seeds[i])
        child_seeds_ints[i] = rng.integers(2**32 - 1)  #change to INT_MAX
    verbose = False    
    reverse = False
    #parallelizing, here the number of workers is set as default to the number of cpus, you can modify it  
    with Pool() as pool:
        res = pool.starmap(nlcd_single, zip(shared_data, [shuffles]*ntrios, [algo]*ntrios, child_seed_ints, [verbose]*ntrios, [reverse]*ntrios))
    
    df=pd.DataFrame(res)
    et = time.time()
    elapsed_time = et - st
    print('start time: ', st,' end time: ',et, 'Execution time:', elapsed_time, 'seconds')
    #seed_list=[]
    #for seed in child_seeds:
    #    rng=np.random.default_rng(seed)
    #    sample_seed=rng.integers(2**32 - 1)
    #    seed_list.append(sample_seed)
    #df['seed']=seed_list
    df.columns=['p_final','p_LassocB','p_LassocA|B','p_AassocB|L','p_LindB|A','OS Test 2','OS Test 4','seed']
    return df

def nlcd_single_for_batch(singletriodata, shuffles, algo, sample_seed=None, verbose=False, reverse=False):
    stopifnot(len(singletriodata)==3)
    ncld_single(L=singletriodata[0], A=singletriodata[1], B=singletriodata[2], shuffles, algo, sample_seed, verbose, reverse)

def nlcd_single(L, A, B, shuffles, algo, sample_seed=None, verbose=True, reverse=False):
    if(sample_seed==None):
        #ss = SeedSequence()
        rng = np.random.default_rng()
        sample_seed = rng.integers(2**32 - 1) #change to INT_MAX
    rng=np.random.default_rng(sample_seed)
    print("Seed set at ",sample_seed)
    sample_seed1 = rng.integers(2**32 - 1) #change to INT_MAX
    sample_seed2 = rng.integers(2**32 - 1) #change to INT_MAX
    np.random.seed(sample_seed1)
    random.seed(sample_seed2)
        
    if reverse==False:
        out=combine_tests(L,A,B,shuffles,algo)
    elif reverse=True:
        out=combine_tests(L,B,A,shuffles,algo)
    else:
        print("Invalid entry for reverse parameter")
    if verbose==True:
        print("The final p value is ",out[0])
        print("Test 1 L assoc B ",out[1])
        print("Test 2 L assoc A | B ",out[2])
        print("Test 3 A assoc B | L ",out[3])
        print("Test 4 L ind B | A ",out[4])
        print("Overlap score from Test 2 ",out[5])
        print("Overlap score from Test 4 ",out[6])
        
    out.append(sample_seed)
    return out


