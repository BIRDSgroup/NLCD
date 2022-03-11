library(MRPC)
library(pcalg)
library(bnlearn)

read_data<-function(path)
{
  dataset<- vector("list", 1000)
  con = file(path, "r")
  for (i in 1:1000)
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
  return (dataset)
}
dataset_0<- read_data("../yeast_residual_data_full_1000_gt_2.txt")#gnd truth 0
dataset_1<- read_data("../yeast_residual_data_full_1000_gt_1.txt") #gnd truth 1
choose_algo <- function(name,genedata)
{
  obj<-NULL
  suffStat<- list(C=cor(genedata),n=nrow(genedata))
  if(name=="MRPC")
    obj<- MRPC(genedata,
               suffStat = suffStat,
               GV = 1,
               FDR = 0.05,
               indepTest = 'gaussCItest',
               labels = colnames(genedata),
               FDRcontrol = 'LOND',
               verbose = FALSE)
  else if(name=="PC")
    obj <- pc(suffStat = suffStat,
              indepTest = gaussCItest,
              alpha = 0.05,
              labels = colnames(genedata),
              verbose = FALSE)
  return (as(obj@graph,"matrix"))
}
perform_theTest<-function(dataset,n,test)
{
  finalresult<-vector("list",n)
  for(i in 1:n)
  {
    temp<-as.data.frame(dataset[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))

    res<-choose_algo(test,temp)
    flat<-as.vector(t(res))
    finalresult[[i]]<-flat
  }
    return (finalresult)
}
specifytest<-function(datanumber,test,n)
{
  result<-NULL
  if(datanumber==1)
    result<-perform_theTest(dataset_1,n,test)
  else if(datanumber==0)
    result<-perform_theTest(dataset_0,n,test)
  finalmatrix<-as.data.frame(do.call(rbind,result))
  colnames(finalmatrix)<-c("L->L","L->A","L->B","A->L","A->A","A->B","B->L","B->A","B->B")
  return (finalmatrix)
}

test="PC"
datanumber=0
tests<-c("MRPC","PC")
datanumbers<-c(0,1)
n=1000
for(i in tests)
{
  for(j in datanumbers)
  {
    finalmatrix<-specifytest(j,i,n)
    labels<-colnames(finalmatrix)
    finalmatrix$model<-apply(finalmatrix,1,function(x) paste0(unlist(labels[which(x %in% 1)]),collapse=","))
    filename=paste(i,"truth",j,"numberofdata",n,".csv",sep="_")
    write.csv(finalmatrix,filename,row.names = FALSE)
    
  }
}


