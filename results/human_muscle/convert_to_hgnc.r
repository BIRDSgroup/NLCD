library(biomaRt)
library(dplyr)
convert_hgnc <-function(data)
{
causalpair<-read.table(data,sep=',',header=TRUE,check.names=TRUE)
causalpair$source<-as.character(causalpair$source)
causalpair$target<-as.character(causalpair$target)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# for source genes 
tempsource<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = causalpair$source, mart = mart,useCache = FALSE) 
causalsource<-merge(x = causalpair, y = tempsource, by.x = "source",by.y= "ensembl_gene_id" , all.x = TRUE)
names(causalsource)[names(causalsource) == 'hgnc_symbol'] <- 'hgnc_source'
causalsource$hgnc_source <- ifelse(is.na(causalsource$hgnc_source), causalsource$source, causalsource$hgnc_source)
causalsource$hgnc_source<- ifelse(causalsource$hgnc_source == '', causalsource$source, causalsource$hgnc_source)
# for target genes 
temptarget<-getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = causalsource$target, mart = mart,useCache = FALSE)  
causalfinal<-merge(x = causalsource, y = temptarget, by.x = "target",by.y= "ensembl_gene_id" , all.x = TRUE)
names(causalfinal)[names(causalfinal) == 'hgnc_symbol'] <- 'hgnc_target'
causalfinal$hgnc_target <- ifelse(is.na(causalfinal$hgnc_target), causalfinal$target, causalfinal$hgnc_target)
causalfinal$hgnc_target<- ifelse(causalfinal$hgnc_target == '', causalfinal$target, causalfinal$hgnc_target)
return(causalfinal)
}

write.csv(convert_hgnc("./results/journal/human_muscle/nlcd/causal.csv"),"./results/journal/human_muscle/nlcd/causal_converted.csv",row.names = FALSE,quote=FALSE)
write.csv(convert_hgnc("./results/journal/human_muscle/nlcd/causal_rev.csv"),"./results/journal/human_muscle/nlcd/causal_rev_converted.csv",row.names=FALSE,quote=FALSE)












