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

def nlcd_batch(shared_data,shuffles,algo,sample_seed=None):
    st = time.time()
    if(sample_seed==None):
        ss = SeedSequence()
    else:
        ss = SeedSequence(sample_seed)
    print("Seed Set at ",ss.entropy)
    
    #unique child seeds for each of the sample 
    child_seeds = ss.spawn(len(shared_data))
    #parallelizing , here the number of workers is set as default to the number of cpus, you can modify it  
    with Pool() as pool:
        res=pool.starmap(main_call, zip(range(len(shared_data)), child_seeds,[shared_data]*len(shared_data),[shuffles]*len(shared_data),[algo]*len(shared_data)))
    
    df=pd.DataFrame(res)
    et = time.time()
    elapsed_time = et - st
    print('start time: ', st,' end time: ',et, 'Execution time:', elapsed_time, 'seconds')
    seed_list=[]
    for seed in child_seeds:
        rng=np.random.default_rng(seed)
        sample_seed=rng.integers(2**32 - 1)
        seed_list.append(sample_seed)
    df['seed']=seed_list
    df.columns=['p_final','p_LassocB','p_LassocA|B','p_AassocB|L','p_LindB|A','OS Test 2','OS Test 4','seed']
    return df

def nlcd_single(L,A,B,shuffles,algo,sample_seed=None):
    if(sample_seed==None):
        ss = SeedSequence()
        rng=np.random.default_rng(ss)
        sample_seed=rng.integers(2**32 - 1)
    print("Seed set at ",sample_seed)
    np.random.seed(sample_seed)
    random.seed(sample_seed)
    out=combine_tests(L,A,B,shuffles,algo)
    print("The final p value is ",out[0])
    print("Test 1 L assoc B ",out[1])
    print("Test 2 L assoc A | B ",out[2])
    print("Test 3 A assoc B | L ",out[3])
    print("Test 4 L ind B | A ",out[4])
    print("Overlap score from Test 2 ",out[5])
    print("Overlap score from Test 4 ",out[6])




