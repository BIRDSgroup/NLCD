
##### findr ##########
# not going to use the below pv value approach 
library(findr)
# nrow=1
findr.lib(rs=47)
# ### example ####
# L <- matrix(rbinom(800, 1, 0.5), nrow = nrow)
# A <- 10 * L+5*(1-L)  + matrix(rnorm(800, 0, 1), nrow = nrow)
# B <- A + matrix(rnorm(800, 0, 1), nrow = nrow)
# 
# findr.pijs_gassist_pv(L,A,B,na=1)
# 
# ###################
# data_types<-c("Linear","Sine","Saw","Indp")
# data_samples<-c("300","500","1000")
# for(dt in data_types)
# {
#   for(ds in data_samples)
#   {
# test_sim<-read_data(paste0("./../data/",dt,ds,".txt"),inputs=100) #100 is the number of trios 
# res_list<-list()
# for(i in 1:100)
# {
# #timeseed=as.integer(Sys.time())
# timeseed<-  sample.int(1000000, 1) #findr seed is integer
# findr.lib(rs=timeseed,nth=1)
# L<-t(matrix(as.integer(unlist(test_sim[[i]][1]))));
# A<-t(matrix(unlist(test_sim[[i]][2])));
# B<-t(matrix(unlist(test_sim[[i]][3])));
# res<-findr.pijs_gassist_pv(L,A,B,na=1)
# res$p_trad<-max(res$p1,res$p2,res$p3)
# res$p_findr<-min(max(res$p2,res$p5),res$p4)
# res$seed<-timeseed
# res_list[[i]]<-res
# }
# res_table<-as.data.frame(do.call(rbind, res_list))
# res_table<-apply(res_table,2,as.character)
# write.csv(res_table,file=paste0("./../results/journal/simulation/findr/",dt,ds,"findr",".csv"),row.names = F,quote = F);
# print(paste0(dt,ds," completed"))
# }
# }
# 
# example<-read.csv("./../results/journal/simulation/findr/Linear300findr.csv",header=T)
# example_indp<-read.csv("./../results/journal/simulation/findr/Indp300findr.csv",header=T)


### combining trios and running findr probability function #####
convert_list_to_matrix <- function(test_sim) {
  L <- lapply(test_sim, function(x) as.integer(unlist(x[1])))
  L_mat <- do.call(rbind, L)
  A <- lapply(test_sim, function(x) unlist(x[2]))
  A_mat <- do.call(rbind, A)
  B <- lapply(test_sim, function(x) unlist(x[3]))
  B_mat <- do.call(rbind, B)
  return(list(L_mat, A_mat, B_mat))
}
convert_data_into_mat<-function(dt,ds,inputs)
{
  test_sim<-read_data(paste0("./../data/",dt,ds,".txt"),inputs)
  result<-convert_list_to_matrix(test_sim)
  return (result)
}

convert_mat_multiple<-function(ds,inputs)
{
  L_mat <- matrix(nrow = 0, ncol = as.integer(ds))
  A_mat <- matrix(nrow = 0, ncol = as.integer(ds))
  B_mat <- matrix(nrow = 0, ncol = as.integer(ds))
  for(dt in c("Sine","Saw","Linear","Indp"))
  {
    result<-convert_data_into_mat(dt,ds,inputs)
    L_mat <-rbind(L_mat, result[[1]])
    A_mat <- rbind(A_mat,result[[2]])
    B_mat <- rbind(B_mat,result[[3]])
  }
  return(list(L_mat, A_mat, B_mat))
}
ds="1000";inputs=100;
result<-convert_mat_multiple(ds,inputs)
L_mat <- result[[1]]
A_mat <- result[[2]]
B_mat <- result[[3]]


results_comb<-findr.pij_gassist(L_mat,A_mat,B_mat,na=1)
diagonals<-diag(results_comb)
causal<-diagonals[1:100]
indp<-diagonals[301:400]


library(PRROC)
groundtruth<- c(rep(1, 100), rep(0, 100))
res_list<-c(causal,indp)
#pr_values <- pr.curve(scores.class0 = diagonals, weights.class0 = groundtruth,curve=T)
pr_values <- pr.curve(scores.class0 = causal,scores.class1=indp,curve=T)
pr_values <- pr.curve(scores.class0 = res_list, weights.class0 = groundtruth,curve=T)
# Plot PR curve
plot(pr_values, main = "Precision-Recall Curve", col = "blue")
legend("bottomright", legend = c("PR Curve"), col = "blue", lty = 1)



####looping it over all the sample sizes and saving the result ######
for(ds in c("300","500","1000"))
{
  result<-convert_mat_multiple(ds,100)
  L_mat <- result[[1]]
  A_mat <- result[[2]]
  B_mat <- result[[3]]
  
  
  results_comb<-findr.pij_gassist(L_mat,A_mat,B_mat,na=1)
  diagonals<-diag(results_comb)
  indp<-diagonals[301:400]
  #write.table(indp,file=paste0("./../results/journal/simulation/findr/","Indp",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  write.table(indp,file=paste0("./../results/journal/simulation/findr_sanity/","Indp",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  causal<-diagonals[1:100]
  #write.table(causal,file=paste0("./../results/journal/simulation/findr/","Sine",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  write.table(causal,file=paste0("./../results/journal/simulation/findr_sanity/","Sine",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  causal<-diagonals[101:200]
  #write.table(causal,file=paste0("./../results/journal/simulation/findr/","Saw",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  write.table(causal,file=paste0("./../results/journal/simulation/findr_sanity/","Saw",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  causal<-diagonals[201:300]
  #write.table(causal,file=paste0("./../results/journal/simulation/findr/","Linear",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
  write.table(causal,file=paste0("./../results/journal/simulation/findr_sanity/","Linear",ds,"findr",".csv"),row.names = F,quote = F,col.names ="p_val")
}


### findr yeast ###### 
#check findr_yeast.ipynb for the preprocessing
#causal data
L_mat<-as.matrix(read.table("./../../../../findr/findrfiles/scripts/findr_yeast/L_mat_findr.csv",sep = ',',row.names = 1,header=T))
A_mat<-as.matrix(read.table("./../../../../findr/findrfiles/scripts/findr_yeast/A_mat_findr.csv",sep = ',',row.names = 1,header=T))
B_mat<-as.matrix(read.table("./../../../../findr/findrfiles/scripts/findr_yeast/B_mat_findr.csv",sep = ',',row.names = 1,header=T))
findr.lib(rs=47)
results_comb<-findr.pij_gassist(L_mat,A_mat,B_mat,na=1,nodiag=T)
#checking if findr is deterministic 
#findr.lib(rs=51)
#results_comb1<-findr.pij_gassist(L_mat,A_mat,B_mat,na=1,nodiag=T)
#sum(results_comb==results_comb1)
write.table(results_comb,"./../../../../findr/findrfiles/scripts/findr_yeast/causal_results.csv",sep = ',',quote=F)
#independent data
L_mat<-as.matrix(read.table("./../../../../findr/findrfiles/scripts/findr_yeast/L_mat_findr_indp.csv",sep = ',',row.names = 1,header=T))
A_mat<-as.matrix(read.table("./../../../../findr/findrfiles/scripts/findr_yeast/A_mat_findr_indp.csv",sep = ',',row.names = 1,header=T))
B_mat<-as.matrix(read.table("./../../../../findr/findrfiles/scripts/findr_yeast/B_mat_findr_indp.csv",sep = ',',row.names = 1,header=T))
#findr.pv()
results_comb<-findr.pij_gassist(L_mat,A_mat,B_mat,na=1,nodiag=T)
#results_comb1<-findr.pij_gassist(L_mat,A_mat,B_mat,na=1,nodiag=T)
write.table(results_comb,"./../../../../findr/findrfiles/scripts/findr_yeast/indp_results.csv",sep = ',',quote=F)




#### findr test bed ###
dt="Saw";ds="1000";inputs=100
result<-convert_data_into_mat(dt,ds,inputs)
L_c <- result[[1]]
A_c <- result[[2]]
B_c <- result[[3]]


dt="Indp";ds="1000";inputs=100
result<-convert_data_into_mat(dt,ds,inputs)
L_i <- result[[1]]
A_i <- result[[2]]
B_i <- result[[3]]

L_comb<-rbind(L_c,L_i)
A_comb<-rbind(A_c,A_i)
B_comb<-rbind(B_c,B_i)

results_comb<-findr.pij_gassist(L_comb,A_comb,B_comb,na=1)
diagonals<-diag(results_comb)
causal<-diagonals[101:200]
indp<-diagonals[301:400]

######################

### MRPC testing ###########
yeast_causal<-read_trios("./../../../../findr/findrfiles/scripts/yeastgt_1_wilko1234_ready.txt",inputs=1234)
yeast_indep<-read_trios("./../../../../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt",inputs=1234)
# going to group based on eQTLs and run MRPC for each eQTL group 
