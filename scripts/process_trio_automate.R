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
system(paste0("cut -f $(paste -s -d ',' ",tissue,"geno_mainIndices.txt) ./allgeno/allgeno.012 > ",tissue,"_geno_mat.txt"))
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
# Set to character() for no covariates
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
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)

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





