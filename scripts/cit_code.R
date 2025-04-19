library(cit)
read_data<-function(path,inputs)
{
  dataset<- vector("list", inputs)
  con = file(path, "r")
  if(grepl("yeast",path)==FALSE & grepl("human",path)==FALSE)
    line = readLines(con, n = 1)  #need this for the newer datasets because there is a seed in every dataset
  #con = file(path, "r")
  for (i in 1:inputs)
  {
    line = readLines(con, n = 1)
    line = readLines(con, n = 1)
    line_edit<-unlist(strsplit(line," "))
    l<-as.vector(as.numeric(line_edit))
    line = readLines(con, n = 1)
    line_edit<-unlist(strsplit(line," "))
    a<-as.vector(as.numeric(line_edit))
    line = readLines(con, n = 1)
    line_edit<-unlist(strsplit(line," "))
    b<-as.vector(as.numeric(line_edit))
    dataset[[i]]<-list(l,a,b)
    
  }
  close(con)
  return (dataset)
}
cit_process<- function(path,inputs,reverse=FALSE,seed=NULL,perms=50,outpath=NULL,normal=FALSE)
{
  if(is.null(seed))
  {seed <- sample(.Machine$integer.max, 1)}
  set.seed(seed)
  dataset<- read_data(path,inputs)
  seed_config<-sample(.Machine$integer.max, inputs)
  p_cit<-numeric(inputs)
  p_cit[1:inputs]=-1
  p_TL<- numeric(inputs)
  p_TL[1:inputs]=-1
  p_TG<- numeric(inputs)
  p_TG[1:inputs]=-1
  p_GL<- numeric(inputs)
  p_GL[1:inputs]=-1
  p_Lind<- numeric(inputs)
  p_Lind[1:inputs]=-1
  for(i in 1:inputs)
  {
        temp<-as.data.frame(dataset[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
        L<-temp[,1]
        A<- temp[,2]
        B<- temp[,3]
        if(normal==TRUE)
        {
          unique_values <- unique(L)
          
          for (value in unique_values) {
            indices <- which(L == value)
            mu <- mean(A[indices])
            sigma <- sd(A[indices])
            A[indices] <- (A[indices]) / sigma
          }
          for (value in unique_values) {
            indices <- which(L == value)
            mu <- mean(B[indices])
            sigma <- sd(B[indices])
            B[indices] <- (B[indices]) / sigma
          }
        }
        if(reverse==FALSE)
        t<-cit.cp(L,A,B,rseed = seed_config[i],n.resampl=perms)
        else if(reverse==TRUE)
        { t<-cit.cp(L,B,A,rseed=seed_config[i],n.resampl = perms)}
        
        p_cit[i]<-t[1]
        p_TL[i]<-t[2]
        p_TG[i]<-t[3]
        p_GL[i]<-t[4]
        p_Lind[i]<-t[5]


    
  }
df <- data.frame(p_cit,p_TL,p_TG,p_GL,p_Lind,childseed=seed_config,Parentseed=seed,stringsAsFactors = FALSE)
if(!is.null(outpath))
{write.csv(df,outpath,row.names=FALSE)}
return(df)
}
####### Simulation with permutation = 100 ######
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Sine500cit100perm.csv')
cit_process("./data/Sine300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Sine300cit100perm.csv')
cit_process("./data/Sine1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Sine1000cit100perm.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Saw500cit100perm.csv')
cit_process("./data/Saw300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Saw300cit100perm.csv')
cit_process("./data/Saw1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Saw1000cit100perm.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linear500cit100perm.csv')
cit_process("./data/Linear300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linear300cit100perm.csv')
cit_process("./data/Linear1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linear1000cit100perm.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Indp500cit100perm.csv')
cit_process("./data/Indp300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Indp300cit100perm.csv')
cit_process("./data/Indp1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Indp1000cit100perm.csv')
cit_process("./data/Para500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Para500cit100perm.csv')
cit_process("./data/Para300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Para300cit100perm.csv')
cit_process("./data/Para1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Para1000cit100perm.csv')
cit_process("./data/Paravar500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravar500cit100perm.csv')
cit_process("./data/Paravar300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravar300cit100perm.csv')
cit_process("./data/Paravar1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravar1000cit100perm.csv')
cit_process("./data/Linearvar500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvar500cit100perm.csv')
cit_process("./data/Linearvar300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvar300cit100perm.csv')
cit_process("./data/Linearvar1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvar1000cit100perm.csv')
####### Simulation with permutation = 500 ######
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Sine500cit500perm.csv')
cit_process("./data/Sine300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Sine300cit500perm.csv')
cit_process("./data/Sine1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Sine1000cit500perm.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Saw500cit500perm.csv')
cit_process("./data/Saw300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Saw300cit500perm.csv')
cit_process("./data/Saw1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Saw1000cit500perm.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Linear500cit500perm.csv')
cit_process("./data/Linear300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Linear300cit500perm.csv')
cit_process("./data/Linear1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Linear1000cit500perm.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Indp500cit500perm.csv')
cit_process("./data/Indp300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Indp300cit500perm.csv')
cit_process("./data/Indp1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Indp1000cit500perm.csv')
cit_process("./data/Para500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Para500cit500perm.csv')
cit_process("./data/Para300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Para300cit500perm.csv')
cit_process("./data/Para1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Para1000cit500perm.csv')
cit_process("./data/Paravar500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Paravar500cit500perm.csv')
cit_process("./data/Paravar300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Paravar300cit500perm.csv')
cit_process("./data/Paravar1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Paravar1000cit500perm.csv')
cit_process("./data/Linearvar500.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Linearvar500cit500perm.csv')
cit_process("./data/Linearvar300.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Linearvar300cit500perm.csv')
cit_process("./data/Linearvar1000.txt",inputs=100,perms = 500,outpath = './results/journal/simulation/cit/Linearvar1000cit500perm.csv')
####### 10 runs Linear permutation = 100 ########## 
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun1.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun2.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun3.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun4.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun5.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun6.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun7.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun8.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun9.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Linear500cit100permrun10.csv')

########### Paravar with normalize, only A normalize, function now changed to normalize B also after this ############
cit_process("./data/Paravar500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravarnorm500cit100perm.csv',normal=TRUE)
cit_process("./data/Paravar300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravarnorm300cit100perm.csv',normal=TRUE)
cit_process("./data/Paravar1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravarnorm1000cit100perm.csv',normal=TRUE)
cit_process("./data/Linearvar500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvarnorm500cit100perm.csv',normal=TRUE)
cit_process("./data/Linearvar300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvarnorm300cit100perm.csv',normal=TRUE)
cit_process("./data/Linearvar1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvarnorm1000cit100perm.csv',normal=TRUE)

########### Paravar with normalize, both normalize############
cit_process("./data/Paravar500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravarnormboth500cit100perm.csv',normal=TRUE)
cit_process("./data/Paravar300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravarnormboth300cit100perm.csv',normal=TRUE)
cit_process("./data/Paravar1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Paravarnormboth1000cit100perm.csv',normal=TRUE)
cit_process("./data/Linearvar500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvarnormboth500cit100perm.csv',normal=TRUE)
cit_process("./data/Linearvar300.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvarnormboth300cit100perm.csv',normal=TRUE)
cit_process("./data/Linearvar1000.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/Linearvarnormboth1000cit100perm.csv',normal=TRUE)



####### 10 runs Linear permutation = 500 ########## 
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun1.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun2.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun3.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun4.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun5.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun6.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun7.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun8.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun9.csv')
cit_process("./data/Linear500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Linear500cit500permrun10.csv')


######### 10 runs Indp permutation = 100 ##########
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun1.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun2.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun3.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun4.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun5.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun6.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun7.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun8.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun9.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Indp500cit100permrun10.csv')

######### 10 runs Indp permutation = 500 ##########
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun1.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun2.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun3.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun4.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun5.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun6.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun7.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun8.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun9.csv')
cit_process("./data/Indp500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Indp500cit500permrun10.csv')

######### 10 runs Sine permutation = 100 ##########
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun1.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun2.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun3.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun4.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun5.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun6.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun7.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun8.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun9.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Sine500cit100permrun10.csv')

######### 10 runs Sine permutation = 500 ##########
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun1.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun2.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun3.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun4.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun5.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun6.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun7.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun8.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun9.csv')
cit_process("./data/Sine500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Sine500cit500permrun10.csv')

######### 10 runs Saw permutation = 100 ##########
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun1.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun2.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun3.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun4.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun5.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun6.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun7.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun8.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun9.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 100,outpath = './results/journal/10run/runvariation/cit/Saw500cit100permrun10.csv')

######### 10 runs Saw permutation = 500 ##########
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun1.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun2.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun3.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun4.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun5.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun6.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun7.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun8.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun9.csv')
cit_process("./data/Saw500.txt",inputs=100,perms = 500,outpath = './results/journal/10run/runvariation/cit/Saw500cit500permrun10.csv')


####### 10 runs Linear data permutation = 100 ########## 
cit_process("./data/10rundata/Linear500run1.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun1.csv')
cit_process("./data/10rundata/Linear500run2.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun2.csv')
cit_process("./data/10rundata/Linear500run3.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun3.csv')
cit_process("./data/10rundata/Linear500run4.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun4.csv')
cit_process("./data/10rundata/Linear500run5.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun5.csv')
cit_process("./data/10rundata/Linear500run6.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun6.csv')
cit_process("./data/10rundata/Linear500run7.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun7.csv')
cit_process("./data/10rundata/Linear500run8.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun8.csv')
cit_process("./data/10rundata/Linear500run9.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun9.csv')
cit_process("./data/10rundata/Linear500run10.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Linear500cit100permrun10.csv')

####### 10 runs Linear data permutation = 500 ########## 
cit_process("./data/10rundata/Linear500run1.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun1.csv')
cit_process("./data/10rundata/Linear500run2.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun2.csv')
cit_process("./data/10rundata/Linear500run3.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun3.csv')
cit_process("./data/10rundata/Linear500run4.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun4.csv')
cit_process("./data/10rundata/Linear500run5.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun5.csv')
cit_process("./data/10rundata/Linear500run6.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun6.csv')
cit_process("./data/10rundata/Linear500run7.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun7.csv')
cit_process("./data/10rundata/Linear500run8.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun8.csv')
cit_process("./data/10rundata/Linear500run9.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun9.csv')
cit_process("./data/10rundata/Linear500run10.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Linear500cit500permrun10.csv')

####### 10 runs Indp data permutation = 100 ########## 
cit_process("./data/10rundata/Indp500run1.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun1.csv')
cit_process("./data/10rundata/Indp500run2.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun2.csv')
cit_process("./data/10rundata/Indp500run3.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun3.csv')
cit_process("./data/10rundata/Indp500run4.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun4.csv')
cit_process("./data/10rundata/Indp500run5.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun5.csv')
cit_process("./data/10rundata/Indp500run6.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun6.csv')
cit_process("./data/10rundata/Indp500run7.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun7.csv')
cit_process("./data/10rundata/Indp500run8.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun8.csv')
cit_process("./data/10rundata/Indp500run9.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun9.csv')
cit_process("./data/10rundata/Indp500run10.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Indp500cit100permrun10.csv')

####### 10 runs Indp data permutation = 500 ########## 
cit_process("./data/10rundata/Indp500run1.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun1.csv')
cit_process("./data/10rundata/Indp500run2.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun2.csv')
cit_process("./data/10rundata/Indp500run3.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun3.csv')
cit_process("./data/10rundata/Indp500run4.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun4.csv')
cit_process("./data/10rundata/Indp500run5.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun5.csv')
cit_process("./data/10rundata/Indp500run6.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun6.csv')
cit_process("./data/10rundata/Indp500run7.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun7.csv')
cit_process("./data/10rundata/Indp500run8.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun8.csv')
cit_process("./data/10rundata/Indp500run9.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun9.csv')
cit_process("./data/10rundata/Indp500run10.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Indp500cit500permrun10.csv')


####### 10 runs Sine data permutation = 100 ########## 
cit_process("./data/10rundata/Sine500run1.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun1.csv')
cit_process("./data/10rundata/Sine500run2.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun2.csv')
cit_process("./data/10rundata/Sine500run3.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun3.csv')
cit_process("./data/10rundata/Sine500run4.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun4.csv')
cit_process("./data/10rundata/Sine500run5.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun5.csv')
cit_process("./data/10rundata/Sine500run6.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun6.csv')
cit_process("./data/10rundata/Sine500run7.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun7.csv')
cit_process("./data/10rundata/Sine500run8.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun8.csv')
cit_process("./data/10rundata/Sine500run9.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun9.csv')
cit_process("./data/10rundata/Sine500run10.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Sine500cit100permrun10.csv')

####### 10 runs Sine data permutation = 500 ########## 
cit_process("./data/10rundata/Sine500run1.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun1.csv')
cit_process("./data/10rundata/Sine500run2.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun2.csv')
cit_process("./data/10rundata/Sine500run3.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun3.csv')
cit_process("./data/10rundata/Sine500run4.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun4.csv')
cit_process("./data/10rundata/Sine500run5.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun5.csv')
cit_process("./data/10rundata/Sine500run6.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun6.csv')
cit_process("./data/10rundata/Sine500run7.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun7.csv')
cit_process("./data/10rundata/Sine500run8.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun8.csv')
cit_process("./data/10rundata/Sine500run9.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun9.csv')
cit_process("./data/10rundata/Sine500run10.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Sine500cit500permrun10.csv')

####### 10 runs Saw data permutation = 100 ########## 
cit_process("./data/10rundata/Saw500run1.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun1.csv')
cit_process("./data/10rundata/Saw500run2.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun2.csv')
cit_process("./data/10rundata/Saw500run3.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun3.csv')
cit_process("./data/10rundata/Saw500run4.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun4.csv')
cit_process("./data/10rundata/Saw500run5.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun5.csv')
cit_process("./data/10rundata/Saw500run6.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun6.csv')
cit_process("./data/10rundata/Saw500run7.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun7.csv')
cit_process("./data/10rundata/Saw500run8.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun8.csv')
cit_process("./data/10rundata/Saw500run9.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun9.csv')
cit_process("./data/10rundata/Saw500run10.txt",inputs=100,perms = 100,outpath = './results/journal/10run/datavariation/cit/Saw500cit100permrun10.csv')

####### 10 runs Saw data permutation = 500 ########## 
cit_process("./data/10rundata/Saw500run1.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun1.csv')
cit_process("./data/10rundata/Saw500run2.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun2.csv')
cit_process("./data/10rundata/Saw500run3.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun3.csv')
cit_process("./data/10rundata/Saw500run4.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun4.csv')
cit_process("./data/10rundata/Saw500run5.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun5.csv')
cit_process("./data/10rundata/Saw500run6.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun6.csv')
cit_process("./data/10rundata/Saw500run7.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun7.csv')
cit_process("./data/10rundata/Saw500run8.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun8.csv')
cit_process("./data/10rundata/Saw500run9.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun9.csv')
cit_process("./data/10rundata/Saw500run10.txt",inputs=100,perms = 500,outpath = './results/journal/10run/datavariation/cit/Saw500cit500permrun10.csv')




##### yeast data fresh, from findr generated ########



### yeast updated data 1234trios perms=100###
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko100/cit/yeast_causal_1234.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko100/cit/yeast_indp_1234.csv',perms=100)

### yeast updated data 1234trios perms=500###
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko500/cit/yeast_causal_1234.csv',perms=500)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko500/cit/yeast_indp_1234.csv',perms=500)

#### 10 runs yeast cit causal perms = 100 ##### 
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_1.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_2.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_3.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_4.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_5.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_6.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_7.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_8.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_9.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_causal_10.csv',perms=100)

#### 10 runs yeast cit indep perms = 100 #### 
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_1.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_2.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_3.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_4.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_5.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_6.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_7.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_8.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_9.csv',perms=100)
cit_process("./../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234,outpath = './results/journal/yeast/wilko10runs/cit/yeast_indp_10.csv',perms=100)


##### ignore the below lines, older analaysis ####

### model specificity ######### 
cit_process("./data/model_specific.txt",inputs=100,perms = 100,outpath = './results/model_specific_cit100.csv')
cit_process("./data/model_specific.txt",inputs=100,perms = 500,outpath = './results/model_specific_cit500.csv')

unique_values <- unique(L)

for (value in unique_values) {
  indices <- which(L == value)
  sigma <- sd(A[indices])
  A[indices] <- (A[indices] ) / sigma
}
for (value in unique_values) {
  indices <- which(L == value)
  sigma <- sd(B[indices])
  B[indices] <- (B[indices] ) / sigma
}
temp=read_data("./data/Para500.txt",100)
tempo<-as.data.frame(temp[[4]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L<-tempo[,1]
A<- tempo[,2]
B<- tempo[,3]

mus<-read_data("./gtex/human_muscle_deseq.txt",inputs=5239)
mus1<-read_data("./gtex/muscle/human_muscle.txt",inputs=5239)
for( i in 1:5239)
{
  if(length(unlist(mus[[i]][1]))==0)
  {
    print(i)
    break
  }
}


temp<-as.data.frame(dataset[[5]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
dataset<-read_data("./data/Paravar500.txt",inputs=10)
