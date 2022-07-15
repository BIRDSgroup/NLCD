library(cit)
#inputs=10000
#inputs=62296
#inputs=209157
inputs=62296
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
count=0
a=c()
for(i in 1:inputs)
{
  print(i)
      temp<-as.data.frame(dataset_yeast[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
      L<-temp[,1]
      A<- temp[,2]
      B<- temp[,3]
      if(var(A[L==0])>1)
      {
        print("in")
        if(cit.cp(L,A,B)[1]<0.05)
        {print("causal")
          a=c(a,i)
          count=count+1;
          }
      }
}
      