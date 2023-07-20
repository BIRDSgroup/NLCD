library(dplyr)

path.tis1.log2.expr.resi <- './gtex/Muscle-Skeletal/tissue.mrna.log2.DESeqnormalized.rediual.rds'
path.tis2.log2.expr.resi <- path.tis1.log2.expr.resi
path.tis1.geno <- './gtex/Muscle-Skeletal/Muscle Skeletal/matrixeqtl/genotype.txt'
path.trios <- './gtex/Muscle-Skeletal/Muscle Skeletal/trios.rds'
path.tiss1.me <- './gtex/Muscle-Skeletal/me.var.log2.rds'
path.tiss2.me <- path.tiss1.me


tis1.log2.expr.resi <- readRDS(path.tis1.log2.expr.resi)
tis2.log2.expr.resi <- readRDS(path.tis2.log2.expr.resi)

path.tis1.geno

genotype <- read.table(path.tis1.geno, header = T, strip.white = T,
                       stringsAsFactors = F, check.names = F)
rownames(genotype) <- genotype$id
genotype <- genotype[, -c(1)]

table(genotype == './.')
genotype[(genotype == './.')] <- NA
table(is.na(genotype))
path.trios
trios <- readRDS(path.trios)



head(tis1.log2.expr.resi[, 1:4])  
head(tis2.log2.expr.resi[, 1:4])
head(trios)
head(genotype[, 1:4])

# align all samples 
indv.tis1 <- colnames(tis1.log2.expr.resi)
indv.tis2 <- colnames(tis2.log2.expr.resi)
indv.geno <- colnames(genotype)
indv.com <- intersect(indv.tis1, intersect(indv.tis2, indv.geno))



length(indv.tis1)


length(indv.tis2)
length(indv.geno)
length(indv.com)

tis1.log2.expr.resi <- tis1.log2.expr.resi[, indv.com]


tis2.log2.expr.resi <- tis2.log2.expr.resi[, indv.com]
genotype <- genotype[, indv.com]

dim(tis1.log2.expr.resi)
dim(tis2.log2.expr.resi)
dim(genotype)

table(trios$gene_id.x %in% rownames(tis1.log2.expr.resi))
table(trios$gene_id.y %in% rownames(tis2.log2.expr.resi))
table(trios$variant_id %in% rownames(genotype))


rownames(trios) <- 1:nrow(trios)
dim(trios)
filename="human_muscle.txt"
for (i in 1:nrow(trios[ , ])) {
  L_name<- trios$variant_id[i]
  A_name<- trios$gene_id.x[i]
  B_name<- trios$gene_id.y[i]
  
  L <- as.integer(genotype[trios$variant_id[i], ])
  indv.not.na <- !is.na(L)
  L <- L[indv.not.na]
  
  A <- as.numeric(tis1.log2.expr.resi[trios$gene_id.x[i], indv.not.na])
  B <- as.numeric(tis2.log2.expr.resi[as.character(as.character(trios$gene_id.y[i])), indv.not.na])
  
  
  write.table(paste(L_name,A_name,B_name), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(rbind(L), file = filename, row.names = FALSE, col.names = FALSE,append=TRUE)
  write.table(rbind(A), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
  write.table(rbind(B), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE)
}

### checking if L has single variance ####
v<-c()
for(i in 1:3656)
{
  print(i)
  v<-c(v,length(unique(unlist(muscle[[i]][1]))))
}
# none of them are only one value ## 
