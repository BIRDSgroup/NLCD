import time 
st = time.time()
from nlcd_main import *
# first argument is for the algorithm
algo=sys.argv[1]
# second argument is for the data file 
data= sys.argv[2] 
# third argument is for the output file 
outputname= sys.argv[3]
# fourth argument is for the number of permutations 
shuffles=int(sys.argv[4])
# fifth argument is whether to check the trios in reverse i.e. L-> B -> A 
reverse=int(sys.argv[5])
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
dataset = [i for i in zip(L,A,B)]
fo.close()

def main_call(i,child_seed):
    '''
    this function is called for each of the samples using Pooling, returns the sample index,p-values and overlap score
    '''
    rng=np.random.default_rng(child_seed)
    sample_seed1=rng.integers(2**32 - 1)
    sample_seed2=rng.integers(2**32 - 1)
    #child_seed = SeedSequence(child_seed_entropy);  seeds = child_seed.spawn(2)
    np.random.seed(seed1)
    random.seed(seed2)
    #set_random_seed(seed3)  #if set_random_seed needed for ANN
    A=np.array(dataset[i][1])
    B=np.array(dataset[i][2])
    L=np.array(dataset[i][0])
    if(reverse==0):
        return combine_tests(L,A,B,shuffles,algo)
    elif(reverse==1):
        return combine_tests(L,B,A,shuffles,algo)
    else:
        print("Invalid entry for reverse parameter")


if __name__ == '__main__':
    #generate a seed and save it to the start of the output file 
    ss = SeedSequence()
    f=open(outputname,"a")
    f.write("Seed "+ str(ss.entropy))
    f.write("\n")
    f.close()
    #unique child seeds for each of the sample 
    child_seeds = ss.spawn(len(dataset))
    #parallelizing , here the number of workers is set as default to the number of cpus, you can modify it  
    with Pool( initargs=(dataset,)) as pool:
        res=pool.starmap(main_call, zip(range(len(dataset)), child_seeds))
    
    df=pd.DataFrame(res)
    
    et = time.time()
    elapsed_time = et - st
    print("Pooling of "+outputname)
    print('start time: ', st,' end time: ',et, 'Execution time:', elapsed_time, 'seconds')
    df.to_csv(outputname, header=False, index=False,mode='a')
    
