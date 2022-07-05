library(cit)
#inputs=10000
#inputs=62296
#inputs=209157
inputs=121
read_data<-function(path)
{
  dataset<- vector("list", inputs)
  con = file(path, "r")
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
dataset_linear<- read_data("../testing_writingvalues_humans.txt")
#dataset_yeast<- read_data("../../../yeast_full_data/yeast_residual_data_full_209k_gt2.txt")
#indices_used<-read_pickle_file("../../thirdApproach/indicesUsedIndependent.pkl")
#adding +1 since R is 1 indexed 
#indices_used<-indices_used+1
#dataset_linear<- dataset_yeast[indices_used[1:10000]]
#inputs=10000
inputs=121
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
p_res<-character(inputs)
p_res[1:inputs]="NA"
for(i in 1:inputs)
{
  print(i)
  try(
  {
temp<-as.data.frame(dataset_linear[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L<-temp[,1]
A<- temp[,2]
B<- temp[,3]

t<-cit.cp(L,A,B)

#t<-cit.cp(L,B,A)
p_cit[i]<-t[1]
p_TL[i]<-t[2]
p_TG[i]<-t[3]
p_GL[i]<-t[4]
p_Lind[i]<-t[5]

if( t[1]<0.05)
{
p_res[i]<-"causal"
}
else
p_res[i]<-" "

  }
  )}
save(dataset_linear,file="yeast10k.Rdata")
df <- data.frame(p_cit,p_TL,p_TG,p_GL,p_Lind,p_res)
write.csv(df,"result_humans.csv",row.names=FALSE)

