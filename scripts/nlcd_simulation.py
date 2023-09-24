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
    parser.add_argument('-n',"--normal",action='store_true',help = " If you want to normalize A|L default false ")
    parser.add_argument("--seed",type=int,help='Seed to reproduce the result')
    args = parser.parse_args()
    dataset=read_data(args.inputpath)
    #df=nlcd_batch(dataset, args.shuffles, args.algo, args.reverse,args.seed,args.normal)
    p_df,t1loss,t2loss,t4loss,t3loss_0,t3loss_1,t3loss_2=nlcd_batch(dataset,args.shuffles,args.algo,args.reverse,args.seed,args.normal)
    p_df.to_csv(args.outputpath, header=True, index=False,float_format='%f')
    t1loss.to_csv("t1loss.csv",header=True,index=False,float_format='%f')
    t2loss.to_csv("t2loss.csv",header=True,index=False,float_format='%f')
    t4loss.to_csv("t4loss.csv",header=True,index=False,float_format='%f')
    t3loss_0.to_csv("t3loss_0.csv",header=True,index=False,float_format='%f')
    t3loss_1.to_csv("t3loss_1.csv",header=True,index=False,float_format='%f')
    t3loss_2.to_csv("t3loss_2.csv",header=True,index=False,float_format='%f')
