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
# change it to yeastgt_0 for the indpendent dataset
dataset_yeast<-read_data('yeastgt_1_wilko1234_ready.txt',inputs=1234)
inputs=1234
mi_causal<- numeric(inputs)
bcmi_causal<- numeric(inputs)
zvalue_causal<- numeric(inputs)
cor_causal<-numeric(inputs)
spear_causal<-numeric(inputs)
for(i in 1:1234)
{
  print(i)
  try(
    {
      temp<-as.data.frame(dataset_yeast[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
      L<-temp[,1]
      A<- temp[,2]
      B<- temp[,3]
      
      mp<-cmi.pw(A,B)
      mi_causal[i]<-mp$mi
      bcmi_causal[i]<-mp$bcmi
      zvalue_causal[i]<-mp$zvalue
      cor_causal[i]<- cor(A,B)
      spear_causal[i]<-cor(A,B,method = "spearman")
      
      
    }
  )}
df <- data.frame(mi_causal,bcmi_causal,zvalue_causal,cor_causal,spear_causal)
write.csv(df,"mpmicorspear_causal_wilko1234.csv",row.names=FALSE)

