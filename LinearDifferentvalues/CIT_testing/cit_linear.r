library(cit)
read_data<-function(path)
{
  dataset<- vector("list", 5)
  con = file(path, "r")
  for (i in 1:5)
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
dataset_linear<- read_data("testing_writingvalues.txt")
temp<-as.data.frame(dataset_linear[[1]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L<-temp[,1]
A<- temp[,2]
B<- temp[,3]
cit.cp(L,A,B)
