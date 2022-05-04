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
n = 10
gt1 <- subset(yeast, groundtruth == 1)[1:n,]
gt2 <- subset(yeast, groundtruth == 0)[1:n,]

filename = 'yeast_2_residual_cit.txt'
filename1 = 'yeast_residual_data.txt'
Lname = c()
Aname = c()
Bname = c()
for (i in 1:n){
  L <- gsub(":",".",gt1[i,1])
  L <- gsub("/",".",L)
  A_ <- gt1[i,2]
  B_ <- gt1[i,3]
  Lname <- c(Lname,L)
  Aname <- c(Aname,A_)
  Bname <- c(Bname,B_)
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
  write.table(paste('A-->B',L,A_,B_), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(paste('LAB, p_cit:',results_1[1],',','p_TassocL:',results_1[2],',','p_TassocGgvnL:',results_1[3],',','p_GassocLgvnT:',results_1[4],',','p_LindTgvnG:',results_1[5]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(paste('LBA, p_cit:',results_2[1],',','p_TassocL:',results_2[2],',','p_TassocGgvnL:',results_2[3],',','p_GassocLgvnT:',results_2[4],',','p_LindTgvnG:',results_2[5]), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table('\n', file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  
  # write.table('A-->B', file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(l, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(a, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
  # write.table(b, file = "yeast_res_2.txt", row.names = FALSE, col.names = FALSE,append = TRUE)
}

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