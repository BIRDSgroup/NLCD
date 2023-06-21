from nlcd_user import *
from argparse import ArgumentParser
def read_data(data):

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

if __name__ == '__main__':
    parser = ArgumentParser()
    # first argument is for the algorithm
    parser.add_argument('-a','--algo',type=str,help='Algorithm (SVR,KRR,ANN) ',required=True)
    parser.add_argument('-i','--inputpath',type=str,help='Path of the input file',required=True)
    parser.add_argument('-o','--outputpath',type=str,help='Path of the output file',required=True)
    parser.add_argument('-s','--shuffles',type=int,help='Number of permutations',default=100)
    parser.add_argument('-r',"--reverse",action='store_true',help = " If you want to run the test in reverse default false ")
    parser.add_argument("--seed",type=int,help='Seed to reproduce the result')
    args = parser.parse_args()
    dataset=read_data(args.inputpath)
    df=nlcd_batch(dataset, args.shuffles, args.algo, args.reverse,args.seed)
    df.to_csv(args.outputpath, header=False, index=False,mode='a',float_format='%f')
    
