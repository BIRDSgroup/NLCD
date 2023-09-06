library(dplyr)
library(data.table)
#muscle done, adipose now
path.eqtl<-'./GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.egenes.txt.gz'
path.expression<-'./GTEx_Analysis_v8_eQTL_expression_matrices/Adipose_Subcutaneous.v8.normalized_expression.bed.gz'
# line 24183 did not have 33 elements while reading eqtl, hence fill = TRUE 
# read the table egenes file 
eqtl<-read.table(path.eqtl,header=TRUE,check.names = FALSE,fill=TRUE)
# print a column using awk
# awk '{print $1,$3}' <namefile>
# run it on terminal 
#one time all geneotype extraction 
#vcftools --gzvcf ../../../../../../../private-data/GTEX/genotype.v8/dataclient_linux_v8genotype/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.vcf.gz --012 --out allgeno
# https://unix.stackexchange.com/questions/245592/print-common-lines-between-two-files-along-with-their-line-numbers-in-those-files
#awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' file1 file2
# cat allgeno.012.pos | wc -l


#cat allgeno.012.pos | wc -l
# 10770860 rows 
#  awk '{print NF}' allgeno.012 | sort -nu | head -n 1
# 10770861 columns
#  awk '{print NF}' allgeno.012 | sort -nu | tail -n 1
# 10770861 columns 
# awk '{print $1}' allgeno.012 
#https://stackoverflow.com/questions/5761212/how-do-i-count-the-number-of-rows-and-columns-in-a-file-using-bash
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

write.table(unique(eqtl_filtered$variant_id),"muscle_genos.txt",row.names = F,col.names=F,quote=F)
# duplicated indices: 2961  5463  8254  8563  9198 10062 11859 12484 why? 
# genotype values repeated for some reason: check these rows, 2961 and 2960 are repeated at 2133406 and 2144505 
# indices in the allgeno.012.pos file, check why are they repeated; chr3 123479874 repeated
#line no: 2133405:chr3    123479874
#line no: 2133406:chr3    123479874
# they have different allele, zcat the lookup table for this chromosome. C_T or C_A at the end. So
# going to extract the full name from WGS_Feature_overlap file, which has MAF >0.01 satisfied
# Use the below command 
# allgeno.012.pos doesnt have the full name of the chromosome,hence there is a confusion
# going to use the following file for full name ( one time command on the terminal )
#zcat WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt.gz | tail -n +2 | cut -f1 > allvariant_fullname.txt

#checking with the main file for the column mapping, this will give us the column numbers of the main allgeno.012 file 
# the following command is tissue specific
# awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' adipose_genos.txt allvariant_fullname.txt > adiposegeno_mapping.txt
system("awk 'FNR==NR{l[$0]=NR; next}; $0 in l{print $0, l[$0], FNR}' muscle_genos.txt allvariant_fullname.txt > musclegeno_mapping.txt ")

# cutting the column and adding +1 to the indices coz allgeno.012 has the first column as indices, so its 0 indexed 
# cut -d ' ' -f3 adiposegeno_mapping.txt | awk '{print $1 + 1}' > adiposegeno_mainIndices.txt

system("cut -d ' ' -f3 musclegeno_mapping.txt | awk '{print $1 + 1}' > musclegeno_mainIndices.txt")
# going to cut the columns from the mainfile adipose indices

# cut -f $(paste -s -d ',' adiposegeno_mainIndices.txt) ./allgeno/allgeno.012 > adipose_geno_mat.txt

system("cut -f $(paste -s -d ',' musclegeno_mainIndices.txt) ./allgeno/allgeno.012 > muscle_geno_mat.txt")
#  paste allgeno/allgeno.012.indv <(cat adipose_geno_mat.txt) > adipose_geno_mat_row.txt

#the below command does not work in system but it works in terminal 
system("paste allgeno/allgeno.012.indv <(cat muscle_geno_mat.txt) > muscle_geno_mat_row.txt")
#  paste -s -d '\t' <(echo "sample/variant" | cut - -d ' ' adiposegeno_mapping.txt -f1) | cat - adipose_geno_mat_row.txt > adipose_geno_mat_final.txt
# this also works in terminal only 
system("paste -s -d '\t' <(echo \"sample/variant\" | cut - -d ' ' musclegeno_mapping.txt -f1) | cat - muscle_geno_mat_row.txt > muscle_geno_mat_final.txt")
# read the genotype
genotype<-read.table('./adipose/adipose_geno_mat_final.txt',header=F,check.names = FALSE)
# transpose the matrix
genotype<-t(genotype)
colnames(genotype) <- genotype[1,]
genotype <- genotype[-1,]
colnames(genotype)[1]<-"id"
row.names(genotype) <- NULL
genotype<-as.data.frame(genotype)
genotype<-genotype[match(unique(eqtl_filtered$variant_id),genotype$id),]
expression <- expression[, -(1:3)]

cov <- read.table('./GTEx_Analysis_v8_eQTL_covariates/Adipose_Subcutaneous.v8.covariates.txt', 
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
write.table(x = expression, file = "expression.txt",
            sep = '\t', quote = F, row.names = F)

write.table(x = genotype, file = "genotype.txt",
            sep = '\t', quote = F, row.names = F)

write.table(x = cov, file = "covariates.txt",
            sep = '\t', quote = F, row.names = F)

write.table(x = geneloc, file = "geneloc.txt",
            sep = '\t', quote = F, row.names = F)

write.table(x = snpsloc, file = "snpsloc.txt",
            sep = '\t', quote = F, row.names = F)

# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = '.';

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

#using ./. as NA 54269 and 4364
nrow(me$cis$eqtls)
nrow(me$trans$eqtls)

## Plot the Q-Q plot of local and distant p-values

plot(me)
saveRDS(me, file = "muscle_sub_cistrans.rds")
cis_merge<- eqtl_filtered %>% dplyr::select(variant_id,gene_id)
trans_merge<- me$trans$eqtls %>% dplyr::select(variant_id=snps,gene_id=gene)

trios=merge(cis_merge,trans_merge,by.x="variant_id",by.y="variant_id")

table(duplicated(trios %>% dplyr::select(variant_id, gene_id.x, gene_id.y)))

saveRDS(trios, file = "muscle_trios.rds")



### extra analysis ### 


cistable<- eqtl_filtered %>% dplyr::select(variant_id,gene_id,gene_chr,gene_start,gene_end,variant_pos)
cistable$tss<-pmin(abs(cistable$variant_pos - cistable$gene_start), abs(cistable$variant_pos -  cistable$gene_end))
table(cistable$tss < 1e6)

trans <- me$trans$eqtls
trans <- trans %>% dplyr::select(variant_id=snps,gene_id= gene, FDR.y = FDR)

head(cistable)
head(trans)

trans <- merge(trans, eqtl %>% dplyr::select(gene_id, gene_chr, gene_start, gene_end,variant_pos), by.x = 'gene_id', by.y = 'gene_id')

trans$tss <- pmin(abs(trans$variant_pos - trans$gene_start), abs(trans$variant_pos -  trans$gene_end))
table(trans$tss < 1e6)

trios <- merge(cistable, trans, by.x = 'variant_id', by.y = 'variant_id')
head(trios) %>% View()
table(duplicated(trios %>% dplyr::select(variant_id, gene_id.x, gene_id.y)))
thres <- 1e-5
table( trios$FDR.y < thres)
table(as.character(trios$gene_id.x) == as.character(trios$gene_id.y))

##### build the tissue.txt file which contains all the trios with values #######

# run the below command to get the hgcn mapping which will be used later on 
#  zcat gene_reads_2017-06-05_v8_muscle_skeletal.gct.gz | tail -n +3 | cut -f2,3 > muscle_ensg_hgcn.txt

tissue <- read.table('adipose/gene_reads_2017-06-05_v8_adipose_subcutaneous.gct.gz', skip = 2,
                     header = T, strip.white = T, stringsAsFactors = F, check.names = F,row.names = 2)
tissue<-tissue[,-c(1,2)]
colnames(tissue)<-paste(tstrsplit(colnames(tissue), '[-]')[[1]],tstrsplit(colnames(tissue), '[-]')[[2]],sep="-")
# make sure that all the samples present in indv.comm are present in the reads matrix
stopifnot (length(setdiff(indv.comm,colnames(tissue)))==0)
tissue<-tissue[,indv.comm]
stopifnot( colnames(tissue) == colnames(expression[,-c(1)]))
# make sure that all the genes present in the expression matrix are present in the reads matrix
stopifnot (length(setdiff(expression$gene_id,rownames(tissue)))==0)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = tissue, colData = as.data.frame(colnames(tissue)), design = ~ 1)
dds <- estimateSizeFactors(dds)
tissue.mrna.count.DESeqnormalized <- counts(dds, normalized=TRUE)
dim(tissue.mrna.count.DESeqnormalized)
tissue.mrna.log2.DESeqnormalized <- log2(tissue.mrna.count.DESeqnormalized + 0.5)
tissue.subset<- tissue.mrna.log2.DESeqnormalized[expression$gene_id,]
#Step 2 co variate adjustment ######### 
trios<-readRDS('./adipose/adipose_trios.rds')
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

filename="human_adipose_deseq.txt"
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




## backup #######
#Step 2 co variate adjustment ######### 
trios<-readRDS('./muscle/muscle_trios.rds')
expression_adj<-expression[,-c(1)]
cov_adj<-cov[,-c(1)]

library(pbmcapply) # parallel processing of any function 
expression_adjusted <- pbmclapply(1:nrow(expression_adj), FUN = function(x) {
  A <- as.numeric(expression_adj[x, ])
  model <- lm(A ~ t(cov_adj))
  
  A.residual <- residuals(model) + model$coefficients["(Intercept)"]
  
  return(A.residual)    
}, mc.cores = 50, ignore.interactive = F)

expression_adjusted<-as.data.frame(do.call(rbind,expression_adjusted))
rownames(expression_adjusted)<-expression$gene_id
colnames(expression_adjusted)<-colnames(expression_adj)

## build the text file #####
genotype_mod<- genotype
table(genotype_mod==-1)
genotype_mod[genotype_mod==-1]<- NA
table(is.na(genotype_mod))
rownames(genotype_mod)<- genotype_mod[,1]
genotype_mod[,1]<-NULL

#check 
stopifnot(colnames(genotype_mod) == colnames(expression_adjusted))

filename="human_muscle_intercept.txt"
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
  
  A <- as.numeric(expression_adjusted[trios$gene_id.x[i], indv.not.na])
  B <- as.numeric(expression_adjusted[trios$gene_id.y[i], indv.not.na])
  
  stopifnot((length(A)==length(B)) && (length(B)==length(L)))
  
  write.table(paste(L_name,A_name,B_name), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE,quote=F)
  write.table(rbind(L), file = filename, row.names = FALSE, col.names = FALSE,append=TRUE,quote = F)
  write.table(rbind(A), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE,quote = F)
  write.table(rbind(B), file = filename, row.names = FALSE, col.names = FALSE,append = TRUE,quote=F)
}


### process the result ####### 

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

read_trios<-function(path,inputs)
{
  L<-c()
  A<-c()
  B<-c()
  con = file(path, "r")
  if(grepl("yeast",path)==FALSE & grepl("human",path)==FALSE)
    line = readLines(con, n = 1)  #need this for the newer datasets because there is a seed in every dataset
  #con = file(path, "r")
  for (i in 1:inputs)
  {
    line = readLines(con, n = 1)
    line_trios<-unlist(strsplit(line," "))
    L<-c(L,line_trios[1])
    A<-c(A,line_trios[2])
    B<-c(B,line_trios[3])
    line = readLines(con, n = 1)
    line = readLines(con, n = 1)
    line = readLines(con, n = 1)
  }
  close(con)
  df<-data.frame(L=L,A=A,B=B)
  return (df)
}


data<- read_data("./muscle/human_muscle_deseq.txt",inputs = 5239)
resultsBtoA<- read.table("./../results/journal/human_muscle/test_muscle_rev.csv",sep=',',header=T,check.names = F)
resultsAtoB<- read.table("./../results/journal/human_muscle/test_muscle.csv",sep=',',header=T,check.names = F)

cross_map<-read.table("./mappability/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz",sep='\t',header=F)
#check for duplication before version numbers are removed 
sum(duplicated(cross_map[,c("V1","V2")]))
#0 duplication 
#check if there is only one dot in the gene names
stopifnot(grepl("^[^.]*\\.[^.]*$",c(cross_map$V1,cross_map$V2)))
#remove the version numbers 
cross_map$V1<-sub("\\..*$","",cross_map$V1)
cross_map$V2<-sub("\\..*$","",cross_map$V2)
#check if any pairs are duplicated 36757 duplicated 
sum(duplicated(cross_map[,c("V1","V2")]))
#there are duplicates after removing the versions 
#get the duplicated data 
#cross_map_dup<-cross_map[duplicated(cross_map[,c("V1","V2")],fromLast=T),]
#get the number of duplicates, max value of cross mappability, min value 
summary_cross_map<-cross_map %>% group_by(V1,V2) %>% summarize(count=n(),max_value = max(V3),min_value=min(V3))
# summary_cross_map has 28,306,345 entries
# total entries 28343102, subtract 36757 , we get 28306345
# check this gene: ENSG00000214717 grep it and see, there is PAR_Y version 
# used this command: zcat hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz  | grep 'ENSG00000182162'| grep 'ENSG00000169084'
# gave the below result 
#ENSG00000169084.13	ENSG00000182162.10	2
#ENSG00000169084.13	ENSG00000182162.10_PAR_Y	1
#ENSG00000169084.13_PAR_Y	ENSG00000182162.10	1
# hence we have three duplicates 
# samme issue with gene mappability 
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

trionames_map<-read_trios("./muscle/human_muscle_deseq.txt",inputs=5239)
trionames_map$A<-sub("\\..*$","",trionames_map$A)
trionames_map$B<-sub("\\..*$","",trionames_map$B)
trionames_map$id=1:nrow(trionames_map)
#can try match command, will keep the order intact  
triosnames_map$mapA = summary_gene_map[match(triosnames_map$A, summary_gene_map$V1),"max_value"]

trionames_Amerge<-merge(trionames_map,summary_gene_map,by.x='A',by.y='V1')
trionames_Amerge<-trionames_Amerge[order(trionames_Amerge$id),]
colnames(trionames_Amerge)<-c("A","L","B","id","A.count","A.max","A.min")
trionames_bothmerge<-merge(trionames_Amerge,summary_gene_map,by.x='B',by.y='V1')
trionames_bothmerge<-trionames_bothmerge[order(trionames_bothmerge$id),]
colnames(trionames_bothmerge)<-c("B","A","L","id","A.count","A.max","A.min","B.count","B.max","B.min")
trionames_dropped<-trionames_bothmerge[,c("id","L","A","B","A.max","B.max")]
# only 358 entries in trionames_cross_merge
# merge in the reverse direction also by creating a duplicate 
# 656 pairs now
#create a duplicate and merge 
summary_cross_map_rev<-summary_cross_map
summary_cross_map_rev[,c("V1","V2")]<-summary_cross_map_rev[,c("V2","V1")]
summary_cross_map_final<-rbind(summary_cross_map,summary_cross_map_rev)
trionames_cross_merge<-merge(trionames_dropped,summary_cross_map_final,by.x=c("A","B"),by.y=c("V1","V2"),all.x=T)
trionames_cross_merge<-trionames_cross_merge[order(trionames_cross_merge$id),]
#sseparating them into two different tables, one with NAs and the other without
trionames_cm_NA<-trionames_cross_merge[is.na(trionames_cross_merge$count),]
trionames_cm<-trionames_cross_merge[!is.na(trionames_cross_merge$count),]

table(trionames_cm_NA$A.max<1,trionames_cm_NA$B.max<1,useNA="ifany")

table(trionames_cm$A.max<1,trionames_cm$B.max<1,useNA="ifany")

sum(is.na(gene_map$V2)) #1237 

sum(is.na(trionames_cm$A.max))    #0
sum(is.na(trionames_cm$B.max))    #0
sum(is.na(trionames_cm_NA$A.max)) #0
sum(is.na(trionames_cm_NA$B.max)) #0


## going to take only those indices where these values are present ####
mappableindices<- trionames_cm$id
# do the below two only if you want to exclude based on mappability 
resultsAtoB<- resultsAtoB[mappableindices,]
resultsBtoA<- resultsBtoA[mappableindices,]

# without p value adjustment 
addmargins(table(resultsAtoB$p_final<=0.05,resultsBtoA$p_final<=0.05))
## losing out a lot of calls after doing p adjust ## 
## we should not adjust B->A, adjust only A to B, adjusting B to A will give more causal calls ##
# we are relaxing the cutoff if we adjust B->A ##
indices<-which(p.adjust(results$p_final,method="none")<=0.05)
addmargins(table(p.adjust(resultsAtoB$p_final,method="BH")<=0.2,resultsBtoA$p_final<=0.05))
#addmargins(table(p.adjust(resultsAtoB$p_final,method="BH")<=0.2,p.adjust(resultsBtoA$p_final,method="BH")<=0.05))
indices<-which(p.adjust(results$p_final,method="none")<=0.05)


# A->B adjusted p value less than 0.2 and B->A greater than 0.05 for cis->trans causal
# Number of causal pairs 
sum((p.adjust(resultsAtoB$p_final,method="BH")<=0.2 & resultsBtoA$p_final>0.05),na.rm=T)  
# this gives 21 pairs after including mappability 
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
write.csv(causalnames_convertB,"./../results/journal/human_muscle/cistranscausal_map.csv",row.names = F,quote=F)
# B->A adjusted p value less than 0.2 and A->B greater than 0.05 for trans -> cis causal

sum((p.adjust(resultsBtoA$p_final,method="BH")<=0.2 & resultsAtoB$p_final>0.05),na.rm=T)
# if zero dont proceed further
causalindices<-which(p.adjust(resultsBtoA$p_final,method="BH")<=0.2 & resultsAtoB$p_final>0.05)

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
write.csv(causalnames_convertB,"./../results/journal/human_muscle/transciscausal.csv",row.names = F,quote=F)




