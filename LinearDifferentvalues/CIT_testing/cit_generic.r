library(cit)
read_data<-function(path)
{
  dataset<- vector("list", 1)
  con = file(path, "r")
  for (i in 1:1)
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
dataset_linear<- read_data("../testing_writingvalues_debugSine.txt")
p_cit<-c()
p_TL<-c()
p_TG<-c()
p_GL<-c()
p_Lind<-c()
p_res<-c()
for(i in 1:1)
{
temp<-as.data.frame(dataset_linear[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L<-temp[,1]
A<- temp[,2]
B<- temp[,3]
t<-cit.cp(L,A,B)
#t<-cit.cp(L,B,A)
p_cit<-c(p_cit,t[1])
p_TL<-c(p_TL,t[2])
p_TG<-c(p_TG,t[3])
p_GL<-c(p_GL,t[4])
p_Lind<-c(p_Lind,t[5])
if(t[1]<0.05)
{
p_res<-c(p_res,"causal")
}
else
p_res<-c(p_res," ")
}
df <- data.frame(p_cit,p_TL,p_TG,p_GL,p_Lind,p_res)
write.csv(df,"debug_cit_Sinewave.csv",row.names = FALSE)

