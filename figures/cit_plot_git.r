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
get_values<-function(L,A,B)
{
  df=data.frame(L,A,B)
  fit=lm(B~A+L)
  input=data.frame(A,L)
  input0_actual=input
  input0_actual[L==1,2]=0
  input1_actual=input
  input1_actual[L==0,2]=1
  B0=predict(fit,input0_actual)
  B1=predict(fit,input1_actual)
  return (list(B0,B1))
}


######### for putting it in manuscript ######################
df<-read_data("./data/Paravar500.txt",inputs=100)
i=5 # indexed 4 in python so 5 in R 
data<-as.data.frame(df[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L=data[,1]
A=data[,2]
B=data[,3]
res=get_values(L,A,B)
B0=res[[1]]
B1=res[[2]]

# the below code is working 
ggplot(data, aes(A, B, color = factor(L))) +
  geom_point() +
  geom_point(aes(A, B0), col = 'red') +
  geom_point(aes(A, B1), col = 'green') + 
  scale_color_manual(values = c("purple", "yellow","red","green")) +theme(legend.position="none")+
  ggtitle("CIT prediction") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("./results/journal/plots/fig3_example_L_cit.svg",plot=last_plot())
ggsave("./results/journal/plots/fig3_example_L_cit.png",plot=last_plot(),dpi=1200)