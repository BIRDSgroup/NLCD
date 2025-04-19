library(cit)
read_data<-function(path,inputs)
{
  # inputs is the number of trios 
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
cit_process<- function(path,inputs,reverse=FALSE,seed=NULL,perms=50,outpath=NULL)
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

# to reproduce the result say Sine data set with 500 samples and 100 permutations, inputs is the number of trios. 
# inputs = 100 for all the simulated datasets 
# for yeast the input is 1234
cit_process("./data/Sine500.txt",inputs=100,perms = 100,outpath = './results/journal/simulation/cit/test.csv',seed=1280338118)



