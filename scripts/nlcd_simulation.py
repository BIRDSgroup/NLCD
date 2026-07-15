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
    parser.add_argument("-l","--losses",action='store_true',help="Losses are stored by default")
    parser.add_argument("--kernel", type=str,   default=None)
    parser.add_argument("--alpha",  type=float, default=None)
    parser.add_argument("--gamma",  type=float, default=None)
    parser.add_argument("--degree", type=int,   default=None)
    # train/test split mode (R2C4: addresses overfitting concern)
    parser.add_argument("--split", action='store_true',
                        help="Use train/test split for fitting MLE/regressors and evaluating loss on held-out test (R2C4)")
    parser.add_argument("--train-frac", type=float, default=0.8,
                        help="Train fraction when --split is set (default 0.8)")

    args = parser.parse_args()
    hyper = {k: v for k, v in dict(kernel=args.kernel, alpha=args.alpha,
                                  gamma=args.gamma, degree=args.degree).items() if v is not None}
    if hyper:
        args.algo = {"name": args.algo, "params": hyper}

    if args.split and args.normal:
        parser.error("--normal cannot be combined with --split (would leak test-set sigma into preprocessing)")
    if args.split and not (0.0 < args.train_frac < 1.0):
        parser.error("--train-frac must be in (0, 1)")

    dataset=read_data(args.inputpath)
    #df=nlcd_batch(dataset, args.shuffles, args.algo, args.reverse,args.seed,args.normal)
    p_df,t1loss,t2loss,t4loss,t3loss_0,t3loss_1,t3loss_2=nlcd_batch(
        dataset, args.shuffles, args.algo, args.reverse, args.seed, args.normal,
        split=args.split, train_frac=args.train_frac,
    )
    # going to change the float_format parameter from %f to string, by default it is string so remove the parameter
    p_df.to_csv(args.outputpath, header=True, index=False)
    if args.losses==True:
        t1loss.to_csv("t1loss.csv",header=True,index=False)
        t2loss.to_csv("t2loss.csv",header=True,index=False)
        t4loss.to_csv("t4loss.csv",header=True,index=False)
        t3loss_0.to_csv("t3loss_0.csv",header=True,index=False)
        t3loss_1.to_csv("t3loss_1.csv",header=True,index=False)
        t3loss_2.to_csv("t3loss_2.csv",header=True,index=False)
