library(mpmi)
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
dataset_yeast<-read_data('/home/aravind/Documents/yeast_analysis/yeast_cov_corrected/yeastgt_0_wilko1752_ready.txt',inputs=1752)
inputs=1752
mi_indep<- numeric(inputs)
bcmi_indep<- numeric(inputs)
zvalue_indep<- numeric(inputs)
cor_indep<-numeric(inputs)
spear_indep<-numeric(inputs)
for(i in 1:1752)
{
  print(i)
  try(
    {
      temp<-as.data.frame(dataset_yeast[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
      L<-temp[,1]
      A<- temp[,2]
      B<- temp[,3]
      
      mp<-cmi.pw(A,B)
      mi_indep[i]<-mp$mi
      bcmi_indep[i]<-mp$bcmi
      zvalue_indep[i]<-mp$zvalue
      cor_indep[i]<- cor(A,B)
      spear_indep[i]<-cor(A,B,method = "spearman")
      
      
    }
  )}
df <- data.frame(mi_indep,bcmi_indep,zvalue_indep,cor_indep,spear_indep)
write.csv(df,"/home/aravind/Documents/yeast_analysis/mpmicorspear_indep_wilko1752.csv",row.names=FALSE)

