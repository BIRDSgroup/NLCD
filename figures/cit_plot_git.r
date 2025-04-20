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
library(ggplot2)
df<-read_data("./data/Paravar500.txt",inputs=100)
i=5 # indexed 4 in python so 5 in R 
data<-as.data.frame(df[[i]],check.rows=FALSE,check.cols=FALSE,col.names=c("L","A","B"))
L=data[,1]
A=data[,2]
B=data[,3]
res=get_values(L,A,B)
B0=res[[1]]
B1=res[[2]]
ggplot(data,aes(A,B,color=factor(L)))+geom_point()+geom_point(aes(A,B0),col='red')+geom_point(aes(A,B1),col='green')+
  ggtitle("CIT prediction")

# save B0 and B1 into a csv 

write.csv(data.frame(B0, B1), file = "B_columns_cit_para_example_L.csv", row.names = FALSE)
# the below code is working 
ggplot(data, aes(A, B, color = factor(L))) +
  geom_point() +
  geom_point(aes(A, B0), col = 'red') +
  geom_point(aes(A, B1), col = 'green') + 
  scale_color_manual(values = c("purple", "yellow","red","green"),name='L') +theme(legend.position="none")+
  #ggtitle("CIT prediction") +
  theme(plot.title = element_text(hjust = 0.5))+theme_bw()
ggsave("./results/journal/plots/fig3_example_L_cit.svg",plot=last_plot())
ggsave("./results/journal/plots/fig3_example_L_cit.png",plot=last_plot(),dpi=1200)
ggsave("./results/journal/plots/fig3_example_L_cit.pdf",plot=last_plot(),dpi=1200)



### muscle scatter plot ### 
load("big_table_muscle.rdat")
default_colors <- scale_color_hue()$palette(2)
# Example indices for highlighted points
# irf1-pmse1 and irf1-parp10
highlight_indices <- c(3255, 2944)

# Split the dataset into highlighted and non-highlighted points
highlight_data <- big_table_tissue[rownames(big_table_tissue) %in% highlight_indices, ]
non_highlight_data <- big_table_tissue[!(rownames(big_table_tissue) %in% highlight_indices), ]

# Create a label column for non-highlighted points
non_highlight_data$label <- ifelse(non_highlight_data$ismappable == 1, "Yes", "No")
highlight_data$label <- ifelse(highlight_data$ismappable == 1, "Yes", "No")  # Match their `Yes`/`No` status

# Plot with separate layers
ggplot() +
  # Non-highlighted points
  geom_point(data = non_highlight_data, 
             aes(x = -log10(Test_3_forward + 1e-12), 
                 y = -log10(Test_4_forward + 1e-12), 
                 color = label), 
             alpha = 0.2) +
  # Highlighted points with black outline
  geom_point(data = highlight_data, 
             aes(x = -log10(Test_3_forward + 1e-12), 
                 y = -log10(Test_4_forward + 1e-12), 
                 fill = label), 
             size = 3, 
             shape = 21, 
             color = "black", # Black outline
             stroke = 1.5,alpha=1,show.legend = FALSE) + # Thickness of outline
  geom_vline(xintercept = -log10(0.1)) +
  labs(
    x = expression(-log[10]("test 3 p-value")("Correlation")),
    y = expression(atop("-log"[10]("test 4 p-value"), "\n (Causation)")),
    color = "Is tested?"
    # To ensure the legend uses the correct mapping
  ) +
  scale_color_manual(values = c("Yes" = "#00BFC4", "No" = "#F8766D"), na.translate = FALSE) +
  scale_fill_manual(values = c("Yes" = "#00BFC4", "No" = "#F8766D"), na.translate = FALSE) + 
  theme_bw()

