rm(list = ls())
library(dplyr)
library(data.table)
tissues=c('Whole_Blood','Skin_Not_Sun_Exposed_Suprapubic','Skin_Sun_Exposed_Lower_leg','Artery_Tibial','Thyroid','Nerve_Tibial','Lung','Esophagus_Mucosa')
for (tissue in tissues)
{
  print(paste("Starting Tissue ",tissue))
path.eqtl<-paste0('./GTEx_Analysis_v8_eQTL/',tissue,'.v8.egenes.txt.gz')
path.expression<-paste0('./GTEx_Analysis_v8_eQTL_expression_matrices/',tissue,'.v8.normalized_expression.bed.gz')
eqtl<-read.table(path.eqtl,header=TRUE,check.names = FALSE,fill=TRUE)

# read the expression matrix
expression <- read.csv2(gzfile(path.expression), sep = '\t', 
                        header = T, check.names = F, strip.white = T, stringsAsFactors = F)
# to check if the genes are exactly aligned and matching
# no further processing is required now if they are matching 
stopifnot (eqtl$gene_id==expression$gene_id)
stopifnot (length(eqtl$gene_id)==length(unique(eqtl$gene_id)))
stopifnot (length(expression$gene_id)==length(unique(expression$gene_id)))
#all are unique genes, take the genes where the qval is less than 0.05
eqtl_filtered<-eqtl %>% filter(qval <= 0.05)
#length(unique(eqtl_filtered$variant_id) repeated snps are present

write.table(unique(eqtl_filtered$variant_id),paste0(tissue,"_genos.txt"),row.names = F,col.names=F,quote=F)

#checking with the main file for the column mapping, this will give us the column numbers of the main allgeno.012 file 
# the following command is tissue specific
# awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' adipose_genos.txt allvariant_fullname.txt > adiposegeno_mapping.txt
system(paste0("awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' ",tissue,"_genos.txt allvariant_fullname.txt > ",tissue,"geno_mapping.txt "))
system(paste0("cut -d ' ' -f3 ",tissue,"geno_mapping.txt | awk '{print $1 + 1}' > ",tissue,"geno_mainIndices.txt"))
system(paste0("split -l 10000 ",tissue,"geno_mainIndices.txt tmp_part_"))
system(paste0("for file in tmp_part_*;do cut -f $(paste -s -d ',' $file) ./allgeno/allgeno.012 > $file.txt ; done"))
system(paste0("paste tmp_part_*.txt > ",tissue,"_geno_mat.txt"))
system(paste0("rm tmp_part_*"))
#removing the below command and replacing it with the previous 4 commands. Arugement list too long error.
#system(paste0("cut -f $(paste -s -d ',' ",tissue,"geno_mainIndices.txt) ./allgeno/allgeno.012 > ",tissue,"_geno_mat.txt"))
system(paste0("echo \"paste allgeno/allgeno.012.indv <(cat ",tissue,"_geno_mat.txt) > ",tissue,"_geno_mat_row.txt\" | bash"))
# found solution: https://stackoverflow.com/questions/32773537/brackets-in-system-command-in-r
system(paste0("echo \"paste -s -d '\t' <(echo 'sample/variant' | cut - -d ' ' ",tissue,"geno_mapping.txt -f1) | cat - ",tissue,"_geno_mat_row.txt > ",tissue,"_geno_mat_final.txt\" | bash"))

system(paste0("mkdir temp",tissue))
system(paste0("mv ",tissue,"* temp",tissue,"/"))
system(paste0("mv temp",tissue," ",tissue,"/"))
# read the genotype
genotype<-read.table(paste0('./',tissue,'/',tissue,'_geno_mat_final.txt'),header=F,check.names = FALSE)
# transpose the matrix
genotype<-t(genotype)
colnames(genotype) <- genotype[1,]
genotype <- genotype[-1,]
colnames(genotype)[1]<-"id"
row.names(genotype) <- NULL
genotype<-as.data.frame(genotype)
genotype<-genotype[match(unique(eqtl_filtered$variant_id),genotype$id),]
expression <- expression[, -(1:3)]

cov <- read.table(paste0('./GTEx_Analysis_v8_eQTL_covariates/',tissue,'.v8.covariates.txt'), 
                  header = T, strip.white = T, stringsAsFactors = F, check.names = F)

# keep rows and columns in same align, remove the id column and check the common samples 
indv.exp <- colnames(expression)[-c(1)]
indv.geno <- colnames(genotype)[-c(1)]
indv.cov <- colnames(cov)[-c(1)]
indv.comm <- intersect(indv.exp, intersect(indv.geno, indv.cov))
#common samples 
length(indv.comm)

expression <- expression %>% dplyr::select(gene_id, all_of(indv.comm))
genotype <- genotype %>% dplyr::select(id, all_of(indv.comm))
cov <- cov %>% dplyr::select(ID, all_of(indv.comm))

# take all the genes from gene location 
geneloc <- eqtl %>% dplyr::select(geneid= gene_id, chr=gene_chr, s1=gene_start, s2=gene_end)
snpsloc <- eqtl_filtered %>% dplyr::select(snp=variant_id , chr= chr, pos= variant_pos)
snpsloc<- unique(snpsloc)

# checking the dimensions 
all((colnames(genotype) == colnames(expression))[-c(1)])
all((colnames(genotype) == colnames(cov))[-c(1)])
all(expression$gene_id==geneloc$geneid)

all(geneloc$gene_id == expression$gene_id)
all(as.character(snpsloc$snp)==as.character(genotype$id))

## all the matrices are ready ##
## save them to a file ## 
write.table(x = expression, file = paste0(tissue,"/expression.txt"),
            sep = '\t', quote = F, row.names = F)

write.table(x = genotype, file = paste0(tissue,"/genotype.txt"),
            sep = '\t', quote = F, row.names = F)

write.table(x = cov, file = paste0(tissue,"/covariates.txt"),
            sep = '\t', quote = F, row.names = F)

write.table(x = geneloc, file = paste0(tissue,"/geneloc.txt"),
            sep = '\t', quote = F, row.names = F)

write.table(x = snpsloc, file = paste0(tissue,"/snpsloc.txt"),
            sep = '\t', quote = F, row.names = F)

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

print("Staring MatrixEQTL")
## Location of the package with the data files.
base.dir = paste0('./',tissue)

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = paste(base.dir, "/genotype.txt", sep="");
snps_location_file_name = paste(base.dir, "/snpsloc.txt", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "/expression.txt", sep="");
gene_location_file_name = paste(base.dir, "/geneloc.txt", sep="");

# Covariates file name
covariates_file_name = paste(base.dir, "/covariates.txt", sep="");

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 2e-2;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "-1"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
#show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
#show(me$trans$eqtls)

nrow(me$cis$eqtls)
nrow(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values
svg(paste0("./",tissue,"/",tissue,"_qqplot.svg"))
plot(me)
dev.off()

saveRDS(me, file = paste0("./",tissue,"/",tissue,"_sub_cistrans.rds"))
cis_merge<- eqtl_filtered %>% dplyr::select(variant_id,gene_id)
trans_merge<- me$trans$eqtls %>% dplyr::select(variant_id=snps,gene_id=gene)

trios=merge(cis_merge,trans_merge,by.x="variant_id",by.y="variant_id")

table(duplicated(trios %>% dplyr::select(variant_id, gene_id.x, gene_id.y)))

saveRDS(trios, file = paste0("./",tissue,"/",tissue,"_trios.rds"))

print("Starting DeSeq")
##### build the tissue.txt file which contains all the trios with values #######

# run the below command to get the hgcn mapping which will be used later on 
#  zcat gene_reads_2017-06-05_v8_muscle_skeletal.gct.gz | tail -n +3 | cut -f2,3 > muscle_ensg_hgcn.txt


system(paste0("wget https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/gene_reads/gene_reads_2017-06-05_v8_",tolower(tissue),".gct.gz -P ","./",tissue,"/"))
system(paste0("zcat ./",tissue,'/',"gene_reads_2017-06-05_v8_",tolower(tissue),".gct.gz | tail -n +3 | cut -f2,3 > ","./",tissue,"/",tissue,"_ensg_hgcn.txt"))

tiss <- read.table(paste0(tissue,'/gene_reads_2017-06-05_v8_',tolower(tissue),'.gct.gz'), skip = 2,
                     header = T, strip.white = T, stringsAsFactors = F, check.names = F,row.names = 2)
tiss<-tiss[,-c(1,2)]
colnames(tiss)<-paste(tstrsplit(colnames(tiss), '[-]')[[1]],tstrsplit(colnames(tiss), '[-]')[[2]],sep="-")
# make sure that all the samples present in indv.comm are present in the reads matrix
stopifnot (length(setdiff(indv.comm,colnames(tiss)))==0)
tiss<-tiss[,indv.comm]
stopifnot( colnames(tiss) == colnames(expression[,-c(1)]))
# make sure that all the genes present in the expression matrix are present in the reads matrix
stopifnot (length(setdiff(expression$gene_id,rownames(tiss)))==0)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = tiss, colData = as.data.frame(colnames(tiss)), design = ~ 1)
dds <- estimateSizeFactors(dds)
tissue.mrna.count.DESeqnormalized <- counts(dds, normalized=TRUE)
dim(tissue.mrna.count.DESeqnormalized)
tissue.mrna.log2.DESeqnormalized <- log2(tissue.mrna.count.DESeqnormalized + 0.5)
tissue.subset<- tissue.mrna.log2.DESeqnormalized[expression$gene_id,]
#Step 2 co variate adjustment ######### 
trios<-readRDS(paste0('./',tissue,'/',tissue,'_trios.rds'))
cov_adj<-cov[,-c(1)]


library(pbmcapply) # parallel processing of any function 
expression_adjusted <- pbmclapply(1:nrow(tissue.subset), FUN = function(x) {
  A <- as.numeric(tissue.subset[x, ])
  model <- lm(A ~ t(cov_adj))
    
  A.residual <- residuals(model) + model$coefficients["(Intercept)"]
  
  return(A.residual)    
}, mc.cores = 50, ignore.interactive = F)

expression_adjusted<-as.data.frame(do.call(rbind,expression_adjusted))
rownames(expression_adjusted)<-expression$gene_id
colnames(expression_adjusted)<-colnames(tissue)

## build the text file #####
genotype_mod<- genotype
table(genotype_mod==-1)
genotype_mod[genotype_mod==-1]<- NA
table(is.na(genotype_mod))
rownames(genotype_mod)<- genotype_mod[,1]
genotype_mod[,1]<-NULL

#check 
stopifnot(colnames(genotype_mod) == colnames(expression_adjusted))

print("Making the trios text file")
filename=paste0("./",tissue,"/","human_",tissue,"_deseq.txt")
# 3127 chr22_23754197_T_A_b38 this is working for muscle
# 3128 is not working chr22_23970400_G_A_b38 for muscle
# 3883 onwards not working for adipose
for( i in 1:nrow(trios))
{
  L_name<- as.character(trios$variant_id[i])
  A_name<- as.character(trios$gene_id.x[i])
  B_name<- as.character(trios$gene_id.y[i])
  
  L <- genotype_mod[L_name, ]
  indv.not.na <- !is.na(L)
  L <- L[indv.not.na]
  
  A <- as.numeric(expression_adjusted[A_name, indv.not.na])
  B <- as.numeric(expression_adjusted[B_name, indv.not.na])
  
  stopifnot((length(A)==length(B)) && (length(B)==length(L)))
  
  write.table(paste(L_name,A_name,B_name), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE,quote=F)
  write.table(rbind(L), file = filename, row.names = FALSE, col.names = FALSE,append=TRUE,quote = F)
  write.table(rbind(A), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE,quote = F)
  write.table(rbind(B), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE,quote=F)
}

print(paste(tissue," Done"))
}


## the below code can be safely ignored


######### Part 2 mappability and p value adjustment ###############   

# fix the number of inputs, can it be automatic? 

#tissues<-c("Whole_Blood","Lung","Thyroid","Esophagus_Mucosa","Artery_Tibial","Nerve_Tibial","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg")
#tissues already declared in Part 1


cross_map<-read.table("./mappability/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz",sep='\t',header=F)
#remove the version numbers 
cross_map$V1<-sub("\\..*$","",cross_map$V1)
cross_map$V2<-sub("\\..*$","",cross_map$V2)
#check if any pairs are duplicated 36757 duplicated 
sum(duplicated(cross_map[,c("V1","V2")]))
#there are duplicates after removing the versions 
summary_cross_map<-cross_map %>% group_by(V1,V2) %>% summarize(count=n(),max_value = max(V3),min_value=min(V3))

#create a duplicate and merge 
summary_cross_map_rev<-summary_cross_map
summary_cross_map_rev[,c("V1","V2")]<-summary_cross_map_rev[,c("V2","V1")]
summary_cross_map_final<-rbind(summary_cross_map,summary_cross_map_rev)

gene_map<-read.table("./mappability/hg38_gene_mappability.txt.gz",sep='\t',header = F)
#check the duplicates without removing version numbers
sum(duplicated(gene_map[,c("V1")]))
# check if all the genes have exactly one dot
stopifnot(grepl("^[^.]*\\.[^.]*$",gene_map$V1))
gene_map$V1<-sub("\\..*$","",gene_map$V1)
#check the duplicates after removing version numbers
sum(duplicated(gene_map[,c("V1")])) # 45 duplicates
summary_gene_map<-gene_map %>% group_by(V1) %>% summarize(count=n(),max_value = max(V2),min_value=min(V2))
# summary_gene_map has 58174 entries
# total gene map entries 58219 entries, subtract 45, we get 58174 


tissues<-c("adipose","muscle","Whole_Blood","Lung","Thyroid","Esophagus_Mucosa","Artery_Tibial","Nerve_Tibial","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg")

trios_data <- c(
  "Artery_Tibial" = 5738,
  "Lung" = 5872,
  "muscle" = 5239,
  "adipose" = 6689,
  "Thyroid" = 8554,
  "Nerve_Tibial" = 7462,
  "Esophagus_Mucosa" = 6209,
  "Whole_Blood" = 4664,
  "Skin_Not_Sun_Exposed_Suprapubic" = 7390,
  "Skin_Sun_Exposed_Lower_leg" = 8734
)
found<-c()
for(tissue in tissues)
{
tissue<-"adipose"
#resultsBtoA<- read.table(paste0("./../results/journal/","human_",tissue,"/","test_",tissue,"_rev.csv"),sep=',',header=T,check.names = F)
#resultsAtoB<- read.table(paste0("./../results/journal/","human_",tissue,"/","test_",tissue,".csv"),sep=',',header=T,check.names = F)
resultsBtoA<- read.table(paste0("./../results/adipose_highresolution/reverse/test_adipose_rev.csv"),sep=',',header=T,check.names = F)
resultsAtoB<- read.table(paste0("./../results/adipose_highresolution/normal/test_adipose.csv"),sep=',',header=T,check.names = F)

### EEPET p values ####
t1pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t1loss_noheader.txt",sep = ''))
t2pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t2loss_noheader.txt",sep = ''))
t3_0pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_0_noheader.txt",sep = ''))
t3_1pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_1_noheader.txt",sep = ''))
t3_2pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_2_noheader.txt",sep = ''))
t4pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t4loss_noheader.txt",sep = ''))
ptable<-data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)
ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,ptable$t3_2,na.rm = T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)
resultsAtoB<-ptable
t1pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t1loss_noheader.txt",sep = ''))
t2pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t2loss_noheader.txt",sep = ''))
t3_0pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t3loss_0_noheader.txt",sep = ''))
t3_1pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t3loss_1_noheader.txt",sep = ''))
t3_2pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t3loss_2_noheader.txt",sep = ''))
t4pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t4loss_noheader.txt",sep = ''))
ptable<-data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)
ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,ptable$t3_2,na.rm = T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)
resultsBtoA<-ptable
# Artery Tibial trios: 5738
# Lung trios: 5872
# Muscle trios: 5239
# Adipose trios: 6689 
# Thyroid trios: 8554
# Nerve Tibial: 7462 
# Esophagus Mucosa: 6209 
# Whole Blood: 4664
# Skin_Not_Sun_Exposed_Suprapubic: 7390
# Skin_Sun_Exposed_Lower_leg: 8734 

trionames_map<-read_trios(paste0("./",tissue,"/human_",tissue,"_deseq.txt"),inputs=trios_data[tissue])
trionames_map$A<-sub("\\..*$","",trionames_map$A)
trionames_map$B<-sub("\\..*$","",trionames_map$B)
trionames_map$id=1:nrow(trionames_map)
#can try match command, will keep the order intact  
#triosnames_map$mapA = summary_gene_map[match(triosnames_map$A, summary_gene_map$V1),"max_value"]

trionames_Amerge<-merge(trionames_map,summary_gene_map,by.x='A',by.y='V1')
trionames_Amerge<-trionames_Amerge[order(trionames_Amerge$id),]
colnames(trionames_Amerge)<-c("A","L","B","id","A.count","A.max","A.min")
trionames_bothmerge<-merge(trionames_Amerge,summary_gene_map,by.x='B',by.y='V1')
trionames_bothmerge<-trionames_bothmerge[order(trionames_bothmerge$id),]
colnames(trionames_bothmerge)<-c("B","A","L","id","A.count","A.max","A.min","B.count","B.max","B.min")
#for gene mappability take the minimum
trionames_dropped<-trionames_bothmerge[,c("id","L","A","B","A.min","B.min")]

trionames_cross_merge<-merge(trionames_dropped,summary_cross_map_final,by.x=c("A","B"),by.y=c("V1","V2"),all.x=T)
trionames_cross_merge<-trionames_cross_merge[order(trionames_cross_merge$id),]

#trionames_map$AtoB<- resultsAtoB[trionames_map$id,'p_final']
#trionames_map$BtoA<- resultsBtoA[trionames_map$id,'p_final']

#if either of the first gene is 1 then cross map is 0
trionames_cross_merge_mod<-trionames_cross_merge
#change it to A.min, for mappability
trionames_cross_merge_mod$max_value[trionames_cross_merge_mod$A.min==1 | trionames_cross_merge_mod$B.min==1]<- 0
#883+978+299+656
trionames_cm_final<-trionames_cross_merge_mod[!is.na(trionames_cross_merge_mod$max_value),]
#merge the p values of A to B and B to A to this table
#trionames_cm_final$AtoB<- resultsAtoB[trionames_cm_final$id,'p_final']
#trionames_cm_final$BtoA<- resultsBtoA[trionames_cm_final$id,'p_final']
#write.table(trionames_cm_final,"trionames_cm_final.csv",col.names = T,row.names = F,quote=F,sep=',')
#changing the cutoff from 0.75 to 0.7,changing it to 0.9
trionames_cm_final_filtered<-trionames_cm_final[trionames_cm_final$A.min>=0.75 & trionames_cm_final$B.min>=0.75 & trionames_cm_final$max_value==0 ,]
trionames_cm_final_filtered[match(index,trionames_cm_final_filtered$id,nomatch=0),]

## going to take only those indices where these values are present ####
## going to to remove L where the count of unique value is less than 10 ###
data<-read_data("./adipose/human_adipose_deseq.txt",inputs = 6689)
counts=sapply(data,function(x) min(table(x[[1]])))
#addmargins(table(counts<10))
counts_id<-which(counts>=10)
mappableindices<- trionames_cm_final_filtered$id
##union of both ###
both<-intersect(mappableindices,counts_id)
# do the below two only if you want to exclude based on mappability 
resultsAtoB<- resultsAtoB[both,]
resultsBtoA<- resultsBtoA[both,]
index<-both[p.adjust(resultsAtoB$p_final,method="BH")<=0.2 & resultsBtoA$p_final>0.05]

#below line first adjusts for p and then check the mappability
sum((p.adjust(resultsAtoB$p_final,method="BH")<=0.2 & resultsBtoA$p_final>0.05)[both],na.rm=T)

# add the number of trios after mappability and cross mappability 
# add plot of test statistics 
if(sum((p.adjust(resultsAtoB$p_final,method="BH")<=0.1 & resultsBtoA$p_final>0.05),na.rm=T)>0)
{
  print(paste0("In Tissue ",tissue," found causal"))
  found<-c(found,tissue)
}else 
{
  print(paste0("In Tissue ",tissue," not found causal"))
}

}



# this gives 15 pairs after including mappability 
causalindices<-which(p.adjust(resultsAtoB$p_final,method="BH")<=0.2 & resultsBtoA$p_final>0.05)
#remove the above or below line once mappability is finalised 
# inputs=5239 for muscle, adipose 6689
trionames<-read_trios("./muscle/human_muscle_deseq.txt",inputs=5239)
causalnames<-trionames[causalindices,]
rownames(causalnames)<-NULL
causalnames$order_column <- seq_len(nrow(causalnames))
ensg_hgcn<-read.table("./muscle/muscle_ensg_hgcn.txt",header=T,row.names=NULL)
causalnames_convertA<- merge(causalnames,ensg_hgcn,by.x='A',by.y='Name',all.x=T)
causalnames_convertA <- causalnames_convertA[order(causalnames_convertA$order_column), ]
causalnames_convertB<- merge(causalnames_convertA,ensg_hgcn,by.x='B',by.y='Name',all.x=T)
causalnames_convertB <- causalnames_convertB[order(causalnames_convertB$order_column), ]
causalnames_convertB <- causalnames_convertB[,c("L","A","Description.x","B","Description.y")]
colnames(causalnames_convertB)<-c("L","A_ensg","A_hgcn","B_ensg","B_hgcn")
rownames(causalnames_convertB)<-NULL
stopifnot (causalnames_convertB$A_ensg == causalnames$A)
stopifnot (causalnames_convertB$B_ensg == causalnames$B)
#file name changed to map that includes mappability conditions 
write.csv(causalnames_convertB,"./../results/journal/human_muscle/cistranscausal_map_filtered.csv",row.names = F,quote=F)
# B->A adjusted p value less than 0.2 and A->B greater than 0.05 for trans -> cis causal


################## analysis with high resolution p value #########

tissue<-"adipose"
# going to do for adipose 
trionames_map<-read_trios(paste0("./",tissue,"/human_",tissue,"_deseq.txt"),inputs=trios_data[tissue])
trionames_map$A<-sub("\\..*$","",trionames_map$A)
trionames_map$B<-sub("\\..*$","",trionames_map$B)
trionames_map$id=1:nrow(trionames_map)
#can try match command, will keep the order intact  
#triosnames_map$mapA = summary_gene_map[match(triosnames_map$A, summary_gene_map$V1),"max_value"]

trionames_Amerge<-merge(trionames_map,summary_gene_map,by.x='A',by.y='V1')
trionames_Amerge<-trionames_Amerge[order(trionames_Amerge$id),]
colnames(trionames_Amerge)<-c("A","L","B","id","A.count","A.max","A.min")
trionames_bothmerge<-merge(trionames_Amerge,summary_gene_map,by.x='B',by.y='V1')
trionames_bothmerge<-trionames_bothmerge[order(trionames_bothmerge$id),]
colnames(trionames_bothmerge)<-c("B","A","L","id","A.count","A.max","A.min","B.count","B.max","B.min")
#for gene mappability take the minimum
trionames_dropped<-trionames_bothmerge[,c("id","L","A","B","A.min","B.min")]

trionames_cross_merge<-merge(trionames_dropped,summary_cross_map_final,by.x=c("A","B"),by.y=c("V1","V2"),all.x=T)
trionames_cross_merge<-trionames_cross_merge[order(trionames_cross_merge$id),]

#trionames_map$AtoB<- resultsAtoB[trionames_map$id,'p_final']
#trionames_map$BtoA<- resultsBtoA[trionames_map$id,'p_final']

#if either of the first gene is 1 then cross map is 0
trionames_cross_merge_mod<-trionames_cross_merge
#change it to A.min, for mappability
trionames_cross_merge_mod$max_value[trionames_cross_merge_mod$A.min==1 | trionames_cross_merge_mod$B.min==1]<- 0
#883+978+299+656
trionames_cm_final<-trionames_cross_merge_mod[!is.na(trionames_cross_merge_mod$max_value),]
#merge the p values of A to B and B to A to this table
#trionames_cm_final$AtoB<- resultsAtoB[trionames_cm_final$id,'p_final']
#trionames_cm_final$BtoA<- resultsBtoA[trionames_cm_final$id,'p_final']
#write.table(trionames_cm_final,"trionames_cm_final.csv",col.names = T,row.names = F,quote=F,sep=',')
#changing the cutoff from 0.75 to 0.7,changing it to 0.9, changed to 0.99, now to 1
trionames_cm_final_filtered<-trionames_cm_final[trionames_cm_final$A.min>=0.99 & trionames_cm_final$B.min>=0.99 & trionames_cm_final$max_value==0 ,]


## going to take only those indices where these values are present ####
mappableindices<- trionames_cm_final_filtered$id

t1pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t1loss_noheader.txt",sep = ''))
t2pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t2loss_noheader.txt",sep = ''))
t3_0pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_0_noheader.txt",sep = ''))
t3_1pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_1_noheader.txt",sep = ''))
t3_2pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_2_noheader.txt",sep = ''))
t4pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t4loss_noheader.txt",sep = ''))

ptable=data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)

nlcd_p<-read.table("./../results/adipose_highresolution/normal/test_adipose.csv",sep = ',',header = T)

#create a column for idx
rownames(ptable)<-1:nrow(ptable)
ptable<- ptable[rownames(ptable) %in% mappableindices,]

ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,ptable$t3_2,na.rm = T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)
ptable <- ptable[!is.na(ptable$p_final), ]
ptable$p.adjust<- p.adjust(ptable$p_final,method="BH") 

comm_samp<-intersect(rownames(nlcd_p[nlcd_p$p_LindB.A>0.9,]),rownames(ptable))
write.table(ptable,"ptable.csv",quote = F)
ptable$idx<- rownames(ptable)

idx<-rownames(ptable)
original_p<- read.table("./../results/adipose_highresolution/normal/test_adipose.csv",sep=',',header=T,check.names = F)
original_subset<-original_p[idx,]

data<- read_data("./adipose/human_adipose_deseq.txt",inputs = 6689)
ggplot(triodata,aes(A,B))+geom_point()+facet_wrap(~L)+stat_smooth()

## checking old values ######## 
t1pval<-t(read.table("./../results/adipose_highresolution/normal/backup/results_leftissue/pvalues_t1loss.txt",sep = ''))
t2pval<-t(read.table("./../results/adipose_highresolution/normal/backup/results_leftissue/pvalues_t2loss.txt",sep = ''))
t3_0pval<-t(read.table("./../results/adipose_highresolution/normal/backup/results_leftissue/pvalues_t3loss_0.txt",sep = ''))
t3_1pval<-t(read.table("./../results/adipose_highresolution/normal/backup/results_leftissue/pvalues_t3loss_1.txt",sep = ''))
t3_2pval<-t(read.table("./../results/adipose_highresolution/normal/backup/results_leftissue/pvalues_t3loss_2.txt",sep = ''))
t4pval<-t(read.table("./../results/adipose_highresolution/normal/backup/results_leftissue/pvalues_t4loss.txt",sep = ''))
ptable_old<-data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)
rownames(ptable_old)<-1:nrow(ptable_old)
ptable_old<- ptable_old[rownames(ptable_old) %in% mappableindices,]

ptable_old$t3<-pmin(ptable_old$t3_0,ptable_old$t3_1,ptable_old$t3_2,na.rm = T)
ptable_old$p_final<-pmax(ptable_old$t1,ptable_old$t2,ptable_old$t3,ptable_old$t4)
ptable_old <- ptable_old[!is.na(ptable_old$p_final), ]
ptable_old$p.adjust<- p.adjust(ptable_old$p_final,method="BH") 
write.table(ptable_old,"ptable_old.csv",quote=F)


### for the reverse results ####### 
t1pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t1loss_noheader.txt.txt",sep = ''))
t2pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t2loss_noheader.txt.txt",sep = ''))
t3_0pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t3loss_0_noheader.txt.txt",sep = ''))
t3_1pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t3loss_1_noheader.txt.txt",sep = ''))
t3_2pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t3loss_2_noheader.txt.txt",sep = ''))
t4pval<-t(read.table("./../results/adipose_highresolution/reverse/results/pvalues_t4loss_noheader.txt.txt",sep = ''))

ptable_rev=data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)
rownames(ptable_rev)<-1:nrow(ptable_rev)
ptable_rev<- ptable_rev[rownames(ptable_rev) %in% mappableindices,]

ptable_rev$t3<-pmin(ptable_rev$t3_0,ptable_rev$t3_1,ptable_rev$t3_2,na.rm = T)
ptable_rev$p_final<-pmax(ptable_rev$t1,ptable_rev$t2,ptable_rev$t3,ptable_rev$t4)
ptable_rev <- ptable_rev[!is.na(ptable_rev$p_final), ]
ptable_rev$p.adjust<- p.adjust(ptable_rev$p_final,method="BH") 


### simulation values run1 and run2 ####### 

t1pval<-t(read.table("./../results/simulation_highresolution/run2/results/pvalues_t1loss_noheader.txt",sep = ''))
t2pval<-t(read.table("./../results/simulation_highresolution/run2/results/pvalues_t2loss_noheader.txt",sep = ''))
t3_0pval<-t(read.table("./../results/simulation_highresolution/run2/results/pvalues_t3loss_0_noheader.txt",sep = ''))
t3_1pval<-t(read.table("./../results/simulation_highresolution/run2/results/pvalues_t3loss_1_noheader.txt",sep = ''))
t4pval<-t(read.table("./../results/simulation_highresolution/run2/results/pvalues_t4loss_noheader.txt",sep = ''))
ptable=data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t4=t4pval)
ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,na.rm = T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)
ptable <- ptable[!is.na(ptable$p_final), ]
ptable$p.adjust<- p.adjust(ptable$p_final,method="BH") 

nlcd_p<- read.table("./../results/simulation_highresolution/run2/Linear300s1000.csv",sep=',',header=T,check.names = F)
plot(ptable$t1,nlcd_p$p_LassocB)
plot(ptable$t2,nlcd_p$`p_LassocA|B`)
plot(ptable$t4,nlcd_p$`p_LindB|A`)
plot(ptable$t3,nlcd_p$`p_AassocB|L`)

#read the loss files of run 2 ##
t3loss_0<-read.table("./../results/simulation_highresolution/run2/t3loss_0_noheader.csv",sep = ',')
t3loss_1<-read.table("./../results/simulation_highresolution/run2/t3loss_1_noheader.csv",sep = ',')
compare_and_count <- function(col) {
  # Take the first element of the column
  first_value <- col[1]
  
  # Count the number of elements less than the first value
  count <- sum(col[-1] <= first_value)
  
  # Divide the count by the total number of remaining rows, includes p value correction 
  result <- (count+1) / (length(col)-1 + 2)
  
  return(result)
}

# Apply the function to each column of the dataframe
t3_0_manp <- apply(t3loss_0, 2, compare_and_count)
t3_1_manp <- apply(t3loss_1,2,compare_and_count)
t3_manp<-pmin(t3_0_manp,t3_1_manp)
sum(abs(t3_manp-nlcd_p$`p_AassocB|L`))
plot(t3_0_manp,ptable$t3_0)
plot(t3_1_manp,ptable$t3_1)
plot(ptable$t3,nlcd_p$`p_AassocB|L`);
plot(t3_manp,nlcd_p$`p_AassocB|L`,log="xy")



### checking the values for adipose ###  
t1pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t1loss_noheader.txt",sep = ''))
t2pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t2loss_noheader.txt",sep = ''))
t3_0pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_0_noheader.txt",sep = ''))
t3_1pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_1_noheader.txt",sep = ''))
t3_2pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t3loss_2_noheader.txt",sep = ''))
t4pval<-t(read.table("./../results/adipose_highresolution/normal/results/pvalues_t4loss_noheader.txt",sep = ''))

t3loss_0<-read.table("./../results/adipose_highresolution/normal/t3loss_0_noheader.csv",sep = ',')
t3loss_1<-read.table("./../results/adipose_highresolution/normal/t3loss_1_noheader.csv",sep = ',')
t3loss_2<-read.table("./../results/adipose_highresolution/normal/t3loss_2_noheader.csv",sep = ',')
t1loss<-read.table("./../results/adipose_highresolution/normal/t1loss_noheader.csv",sep = ',')
ptable=data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval, t4=t4pval)
ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,ptable$t3_2,na.rm=T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)

t3_0_manp <- apply(t3loss_0, 2, compare_and_count)
t3_1_manp <- apply(t3loss_1,2,compare_and_count)
t3_2_manp <- apply(t3loss_2, 2, compare_and_count)
t3_manp<-pmin(t3_0_manp,t3_1_manp,t3_2_manp,na.rm = T)
t1_manp<-apply(t1loss, 2, compare_and_count)
nlcd_p<- read.table("./../results/adipose_highresolution/normal/test_adipose.csv",sep=',',header=T,check.names = F)
plot(t3_manp,nlcd_p$`p_AassocB|L`) # scattered plot
plot(ptable$t3_0,t3_0_manp)
plot(ptable$t3_1,t3_1_manp)
plot(ptable$t3_2,t3_2_manp)
idx<-which(abs(t3_manp-nlcd_p$`p_AassocB|L`)>0.001)
#trio number 145 debug 
data<-read_data("./adipose/human_adipose_deseq.txt",inputs = 6689)


#### faulty adipose normal results ####
t1pvale<-t(read.table("./../results/adipose_highresolution/backup/normal/results/pvalues_t1loss_noheader.txt",sep = ''))
t2pvale<-t(read.table("./../results/adipose_highresolution/backup/normal/results/pvalues_t2loss_noheader.txt",sep = ''))
t3_0pvale<-t(read.table("./../results/adipose_highresolution/backup/normal/results/pvalues_t3loss_0_noheader.txt",sep = ''))
t3_1pvale<-t(read.table("./../results/adipose_highresolution/backup/normal/results/pvalues_t3loss_1_noheader.txt",sep = ''))
t3_2pvale<-t(read.table("./../results/adipose_highresolution/backup/normal/results/pvalues_t3loss_2_noheader.txt",sep = ''))
t4pvale<-t(read.table("./../results/adipose_highresolution/backup/normal/results/pvalues_t4loss_noheader.txt",sep = ''))

t3loss_0e<-read.table("./../results/adipose_highresolution/backup/normal/t3loss_0_noheader.csv",sep = ',')
t3loss_1e<-read.table("./../results/adipose_highresolution/backup/normal/t3loss_1_noheader.csv",sep = ',')
t3loss_2e<-read.table("./../results/adipose_highresolution/backup/normal/t3loss_2_noheader.csv",sep = ',')
t1losse<-read.table("./../results/adipose_highresolution/backup/normal/t1loss_noheader.csv",sep = ',')
ptablee=data.frame(t1=t1pvale,t2=t2pvale,t3_0=t3_0pvale,t3_1=t3_1pvale,t3_2=t3_2pvale, t4=t4pvale)
ptablee$t3<-pmin(ptablee$t3_0,ptablee$t3_1,ptablee$t3_2,na.rm=T)
ptablee$p_final<-pmax(ptablee$t1,ptablee$t2,ptablee$t3,ptablee$t4)

t3_0_manpe <- apply(t3loss_0e, 2, compare_and_count)
t3_1_manpe <- apply(t3loss_1e,2,compare_and_count)
t3_2_manpe <- apply(t3loss_2e, 2, compare_and_count)
t3_manpe<-pmin(t3_0_manpe,t3_1_manpe,t3_2_manpe,na.rm = T)
t1_manpe<-apply(t1loss, 2, compare_and_count)
nlcd_pe<- read.table("./../results/adipose_highresolution/backup/normal/test_adipose.csv",sep=',',header=T,check.names = F)
plot(t3_manpe,nlcd_pe$`p_AassocB|L`) # scattered plot
plot(ptablee$t3_0,t3_0_manpe)
idx<-which(abs(t3_manpe-nlcd_pe$`p_AassocB|L`)>0.001)

#### simulation result cross check LinearKRR500s100perm ######
sim_test<-read.table("./../results/simulation_crosscheck/LinearKRR500s100perm_test.txt",sep=',',header=T)
sim<-read.table("./../results/journal/simulation/nlcd/LinearKRR500s100perm.csv",sep=',',header=T)
max(abs(sim_test$p_final-sim$p_final))

### going to remove trios where the number of 
hist(sapply(data,function(x) min(table(x[[1]]))))
counts=sapply(data,function(x) min(table(x[[1]])))
addmargins(table(counts<10))
which(counts<10)


