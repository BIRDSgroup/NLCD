from nlcd_user import *
from argparse import ArgumentParser


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
    df.to_csv(args.outputpath, header=True, index=False,mode='a',float_format='%f')
    
