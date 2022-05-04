library(cit)
library(stringr)
library(readxl)

read_exp <- read.delim('SI_Data_01_expressionValues.txt', header = TRUE, sep = "\t")
read_gen <- read.delim('SI_Data_03_genotypes.txt', header = TRUE, sep = "\t")
first_col <- row.names(read_exp)
genes <- c()
for (i in 1:length(first_col)) {
  genes <- c(genes,substr(first_col[i],1,6))
}
row.names(read_exp) <- genes

yeast <- read.csv('neeraj.trios.csv')
cov <- read_excel('SI_Data_02_covariates.xlsx')

#Check segregants order
col <- cov$segregant
seg <- c()
for (i in 1:length(col)) {
  seg <- c(seg,substr(col[i],1,6))
}
all(seg == row.names(read_exp))

set.seed(0) #42
cov$batch = as.factor(cov$batch)
rows <- sample(nrow(yeast))
yeast <- yeast[rows, ]
n = 1000
gt1 <- subset(yeast, groundtruth == 1)[1:n,]
gt2 <- subset(yeast, groundtruth == 0)[1:n,]
gt1<- yeast[yeast$groundtruth==1,]
gt2<- yeast[yeast$groundtruth==0,]
#write.csv(gt1,file="groundtruth1.csv",row.names=FALSE)
#write.csv(gt2,file="groundtruth2.csv",row.names=FALSE)
#filename = 'yeast_2_residual_cit.txt'
#filename1 = 'yeast_residual_data_full_1000_1.txt'
#Lname = c()
#Aname = c()
#Bname = c()
#for groundtruth 1
p_cit_lab<-vector(mode="numeric",length=n)
p_TassocL_lab<-vector(mode="numeric",length=n)
p_TassocGgvnL_lab<-vector(mode="numeric",length=n)
p_GassocLgvnT_lab<-vector(mode="numeric",length=n)
p_LindTgvnG_lab<-vector(mode="numeric",length=n)
p_cit_lba<-vector(mode="numeric",length=n)
p_TassocL_lba<-vector(mode="numeric",length=n)
p_TassocGgvnL_lba<-vector(mode="numeric",length=n)
p_GassocLgvnT_lba<-vector(mode="numeric",length=n)
p_LindTgvnG_lba<-vector(mode="numeric",length=n)
for (i in 1:n){
  print(i)
  L <- gsub(":",".",gt1[i,1])
  L <- gsub("/",".",L)
  A_ <- gt1[i,2]
  B_ <- gsub("-",".",gt1[i,3])
  
  #Lname <- c(Lname,L)
  #Aname <- c(Aname,A_)
  #Bname <- c(Bname,B_)
  l <- read_gen[L]
  a <- read_exp[A_]

  b <- read_exp[B_]
  
  # for (i in 1:1012){
  #   if(l[i,1] < 0){
  #     l[i,1] = 0
  #   }
  # }
  l <- unlist(l)
  a <- unlist(a)
  #a<- (a-mean(a))/sd(a)
  b <- unlist(b)
  #b<- (b-mean(b))/sd(b)
  l[l==-1]<-0
  A <- data.frame(gene=a, batch=cov$batch, OD_covariate=cov$OD_covariate)
  model_A <- lm(A$gene ~ A$batch + A$OD_covariate, data = A)
  B <- data.frame(gene=b,batch=cov$batch,OD_covariate=cov$OD_covariate)
  model_B <- lm(B$gene ~ B$batch + B$OD_covariate, data = B)
  a = model_A$residuals
  b = model_B$residuals
  ########for cit ######################
  #lab[i] <- cit.cp(l, a, b)[[1]]
  #lba[i] <- cit.cp(l, b, a)[[1]]
  #lab[i] <- cit.cp(l, a, b)
  #lba[i] <- cit.cp(l, b, a)
  result_1 <- cit.cp(l, a, b)
  result_2 <- cit.cp(l, b, a)
  p_cit_lab[i]<-result_1[[1]]
  p_TassocL_lab[i]<-result_1[[2]]
  p_TassocGgvnL_lab[i]<-result_1[[3]]
  p_GassocLgvnT_lab[i]<-result_1[[4]]
  p_LindTgvnG_lab[i]<-result_1[[5]]
  
  p_cit_lba[i]<-result_2[[1]]
  p_TassocL_lba[i]<-result_2[[2]]
  p_TassocGgvnL_lba[i]<-result_2[[3]]
  p_GassocLgvnT_lba[i]<-result_2[[4]]
  p_LindTgvnG_lba[i]<-result_2[[5]]
  
  #print(result_1)
  #print(result_2)
  #print('\n')
  #write cit results
  #write.table(paste('A-->B',L,A_,B_), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('LAB, p_cit:',results_1[1]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('LBA, p_cit:',results_2[1]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table('\n', file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('A-->B',L,A_,B_), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(rbind(l), file = filename1, row.names = FALSE, col.names = FALSE,append=TRUE)
  #write.table(rbind(a), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(rbind(b), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table('A-->B', file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(l, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(a, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(b, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
}
dat<- data.frame(v1=p_cit_lab,v2=p_TassocL_lab,v3=p_TassocGgvnL_lab,v4=p_GassocLgvnT_lab,v5=p_LindTgvnG_lab,v6=p_cit_lba,v7=p_TassocL_lba,v8=p_TassocGgvnL_lba,v9=p_GassocLgvnT_lba,v10=p_LindTgvnG_lba)
colnames(dat)<- c("p_cit_lab","p_TassocL_lab","p_TassocGgvnL_lab","p_GassocLgvnT_lab","p_LindTgvnG_lab","p_cit_lba","p_TassocL_lba","p_TassocGgvnL_lba","p_GassocLgvnT_lba","p_LindTgvnG_lba")
write.table(dat,file="groundtruth1_cit.csv",row.names = FALSE)
######for cit groundtruth0##############
p_cit_lab<-vector(mode="numeric",length=n)
p_TassocL_lab<-vector(mode="numeric",length=n)
p_TassocGgvnL_lab<-vector(mode="numeric",length=n)
p_GassocLgvnT_lab<-vector(mode="numeric",length=n)
p_LindTgvnG_lab<-vector(mode="numeric",length=n)
p_cit_lba<-vector(mode="numeric",length=n)
p_TassocL_lba<-vector(mode="numeric",length=n)
p_TassocGgvnL_lba<-vector(mode="numeric",length=n)
p_GassocLgvnT_lba<-vector(mode="numeric",length=n)
p_LindTgvnG_lba<-vector(mode="numeric",length=n)
for (i in 1:n){
  print(i)
  L <- gsub(":",".",gt2[i,1])
  L <- gsub("/",".",L)
  A_ <- gt2[i,2]
  B_ <- gsub("-",".",gt2[i,3])
  
  #Lname <- c(Lname,L)
  #Aname <- c(Aname,A_)
  #Bname <- c(Bname,B_)
  l <- read_gen[L]
  a <- read_exp[A_]
  
  b <- read_exp[B_]
  
  # for (i in 1:1012){
  #   if(l[i,1] < 0){
  #     l[i,1] = 0
  #   }
  # }
  l <- unlist(l)
  a <- unlist(a)
  #a<- (a-mean(a))/sd(a)
  b <- unlist(b)
  #b<- (b-mean(b))/sd(b)
  l[l==-1]<-0
  A <- data.frame(gene=a, batch=cov$batch, OD_covariate=cov$OD_covariate)
  model_A <- lm(A$gene ~ A$batch + A$OD_covariate, data = A)
  B <- data.frame(gene=b,batch=cov$batch,OD_covariate=cov$OD_covariate)
  model_B <- lm(B$gene ~ B$batch + B$OD_covariate, data = B)
  a = model_A$residuals
  b = model_B$residuals
  ########for cit ######################
  #lab[i] <- cit.cp(l, a, b)[[1]]
  #lba[i] <- cit.cp(l, b, a)[[1]]
  #lab[i] <- cit.cp(l, a, b)
  #lba[i] <- cit.cp(l, b, a)
  result_1 <- cit.cp(l, a, b)
  result_2 <- cit.cp(l, b, a)
  p_cit_lab[i]<-result_1[[1]]
  p_TassocL_lab[i]<-result_1[[2]]
  p_TassocGgvnL_lab[i]<-result_1[[3]]
  p_GassocLgvnT_lab[i]<-result_1[[4]]
  p_LindTgvnG_lab[i]<-result_1[[5]]
  
  p_cit_lba[i]<-result_2[[1]]
  p_TassocL_lba[i]<-result_2[[2]]
  p_TassocGgvnL_lba[i]<-result_2[[3]]
  p_GassocLgvnT_lba[i]<-result_2[[4]]
  p_LindTgvnG_lba[i]<-result_2[[5]]
  
  #print(result_1)
  #print(result_2)
  #print('\n')
  #write cit results
  #write.table(paste('A-->B',L,A_,B_), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('LAB, p_cit:',results_1[1]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('LBA, p_cit:',results_2[1]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table('\n', file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('A-->B',L,A_,B_), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(rbind(l), file = filename1, row.names = FALSE, col.names = FALSE,append=TRUE)
  #write.table(rbind(a), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(rbind(b), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table('A-->B', file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(l, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(a, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(b, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
}
dat<- data.frame(v1=p_cit_lab,v2=p_TassocL_lab,v3=p_TassocGgvnL_lab,v4=p_GassocLgvnT_lab,v5=p_LindTgvnG_lab,v6=p_cit_lba,v7=p_TassocL_lba,v8=p_TassocGgvnL_lba,v9=p_GassocLgvnT_lba,v10=p_LindTgvnG_lba)
colnames(dat)<- c("p_cit_lab","p_TassocL_lab","p_TassocGgvnL_lab","p_GassocLgvnT_lab","p_LindTgvnG_lab","p_cit_lba","p_TassocL_lba","p_TassocGgvnL_lba","p_GassocLgvnT_lba","p_LindTgvnG_lba")
write.table(dat,file="groundtruth0_cit.csv",row.names = FALSE)
#for groundtruth 0
filename1 = 'yeast_residual_data_full_1000_gt_2_norm.txt'
#lab<-vector(mode="numeric",length=n)
#lba<-vector(mode="numeric",length=n)
for (i in 1:n){
  print(i)
  L <- gsub(":",".",gt2[i,1])
  L <- gsub("/",".",L)
  A_ <- gt2[i,2]
  B_ <- gsub("-",".",gt2[i,3])
  
  #Lname <- c(Lname,L)
  #Aname <- c(Aname,A_)
  #Bname <- c(Bname,B_)
  l <- read_gen[L]
  a <- read_exp[A_]
  
  b <- read_exp[B_]
  
  # for (i in 1:1012){
  #   if(l[i,1] < 0){
  #     l[i,1] = 0
  #   }
  # }
  l <- unlist(l)
  a <- unlist(a)
  a<- (a-mean(a))/sd(a)
  b <- unlist(b)
  b<- (b-mean(b))/sd(b)
  l[l==-1]<-0
  A <- data.frame(gene=a, batch=cov$batch, OD_covariate=cov$OD_covariate)
  model_A <- lm(A$gene ~ A$batch + A$OD_covariate, data = A)
  B <- data.frame(gene=b,batch=cov$batch,OD_covariate=cov$OD_covariate)
  model_B <- lm(B$gene ~ B$batch + B$OD_covariate, data = B)
  a = model_A$residuals
  b = model_B$residuals
  results_1 = cit.cp(l, a, b)
  results_2 = cit.cp(l, b, a)
  #lab[i] <- cit.cp(l, a, b)[[1]]
  #lba[i] <- cit.cp(l, b, a)[[1]]
  print(results_1)
  print(results_2)
  #print('\n')
  #write cit results
  #write.table(paste('A-->B',L,A_,B_), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('LAB, p_cit:',results_1[1]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table(paste('LBA, p_cit:',results_2[1]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  #write.table('\n', file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(paste('A--B',L,A_,B_), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(rbind(l), file = filename1, row.names = FALSE, col.names = FALSE,append=TRUE)
  write.table(rbind(a), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(rbind(b), file = filename1, row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table('A-->B', file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(l, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(a, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(b, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
}
dat<- data.frame(v1=lab,v2=lba)
colnames(dat)<- c("lab","lba")
write.table(dat,file="groundtruth0_cit.csv",row.names = FALSE)
for (i in 1:n){
  L <- gsub(":",".",gt2[i,1])
  L <- gsub("/",".",L)
  A_ <- gt2[i,2]
  B_ <- gt2[i,3]
  l <- read_gen[L]
  a <- read_exp[A_]
  b <- read_exp[B_]
  # for (i in 1:1012){
  #   if(l[i,1] < 0){
  #     l[i,1] = 0
  #   }
  # }
  l <- unlist(l)
  a <- unlist(a)
  b <- unlist(b)
  A <- data.frame(gene=a, batch=cov$batch, OD_covariate=cov$OD_covariate)
  model_A <- lm(A$gene ~ A$batch + A$OD_covariate, data = A)
  B <- data.frame(gene=b,batch=cov$batch,OD_covariate=cov$OD_covariate)
  model_B <- lm(B$gene ~ B$batch + B$OD_covariate, data = B)
  a = model_A$residuals
  b = model_B$residuals
  results_1 = cit.cp(l, a, b)
  results_2 = cit.cp(l, b, a)
  print(results_1)
  print(results_2)
  print('\n')
  write.table(paste('A---B',L,A_,B_), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(paste('LAB, p_cit:',results_1[1],',','p_TassocL:',results_1[2],',','p_TassocGgvnL:',results_1[3],',','p_GassocLgvnT:',results_1[4],',','p_LindTgvnG:',results_1[5]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(paste('LBA, p_cit:',results_2[1],',','p_TassocL:',results_2[2],',','p_TassocGgvnL:',results_2[3],',','p_GassocLgvnT:',results_2[4],',','p_LindTgvnG:',results_2[5]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table('\n', file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table('A---B', file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(l, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(a, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(b, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
}

i = 2
L <- gsub(":",".",gt1[i,1])
L <- gsub("/",".",L)
A_ <- gt1[i,2]
B_ <- gt1[i,3]
l <- read_gen[L]
a <- read_exp[A_]
b <- read_exp[B_]
l <- unlist(l)
a <- unlist(a)
b <- unlist(b)
A <- data.frame(gene=a, batch=cov$batch, OD_covariate=cov$OD_covariate)
model_A <- lm(A$gene ~ A$batch + A$OD_covariate, data = A)
B <- data.frame(gene=b, batch=cov$batch, OD_covariate=cov$OD_covariate)
model_B <- lm(B$gene ~ B$batch + B$OD_covariate, data = B)
a = model_A$residuals
b = model_B$residuals
results_1 = cit.cp(l, a, b)
results_2 = cit.cp(l, b, a)
print(results_1)
print(results_2)
print('\n')