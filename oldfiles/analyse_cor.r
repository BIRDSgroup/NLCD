GT1<-read.csv("Gt1_full_corr.csv")
hist(GT1$nlcor)
Gt<-GT1[which(GT1$nlcor>=0.5 & GT1$pearson<0.3 & GT1$pearson >-0.3 ),]
Gt
read_data<-function(path)
{
  dataset<- vector("list", 62296)
  con = file(path, "r")
  for (i in 1:62296)
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
dataset_linear<- read_data("../yeast_full_data/yeast_residual_data_full_62k_gt1.txt")
temp<-as.data.frame(dataset_linear[[2733]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L<-temp[,1]
A<- temp[,2]
B<- temp[,3]