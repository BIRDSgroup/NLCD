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
#tissue names mapped with the actual tissue name. The naming convention was different
#for some folders. Hence this mapping was required 
tissue_name <- c(
  "Artery_Tibial" = "artery",
  "Lung" = "lung",
  "muscle" = "muscle",
  "adipose" = "adipose",
  "Thyroid" = "thyroid",
  "Nerve_Tibial" = "nerve",
  "Esophagus_Mucosa" = "esophagus",
  "Whole_Blood" = "blood",
  "Skin_Not_Sun_Exposed_Suprapubic" = "skinsupra",
  "Skin_Sun_Exposed_Lower_leg" = "skinleg"
)
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

#resultsAtoB<-read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/test_",tissue_name[tissue],".csv"),sep =',',header=T)
#resultsBtoA<-read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/test_",tissue_name[tissue],"_rev.csv"),sep =',',header=T)
#colnames(resultsAtoB)<-c("p_final","t1","t2","t3","t4","oscore","cseed","pseed")
#colnames(resultsBtoA)<-c("p_final","t1","t2","t3","t4","oscore","cseed","pseed")
#save(file="results.muscle.normal.rdat",resultsAtoB,resultsBtoA)
#the below is the epept values, see checkforcausal.R for the commands
#save(file="results.muscle.epept.rdat",resultsAtoB,resultsBtoA) 
finaltable<- function(filterdf,resultsAtoB,resultsBtoA,tissue)
{
  ##union of both ###
  bootresults<-resultsAtoB
  stopifnot(nrow(filterdf)==nrow(resultsAtoB))
  stopifnot(nrow(filterdf)==nrow(resultsBtoA))
  bootresults<-cbind(bootresults,filterdf)
  bootresults$ismappable<-bootresults$mapA>=0.75 & bootresults$mapB>=0.75 & bootresults$cmap==0 
  bootresults$ismappablewcount<-bootresults$mapA>=0.75 & bootresults$mapB>=0.75 & bootresults$cmap==0 & bootresults$counts>=10 
  #converting the gene name form ensg format to hgcn format
  ensg_hgcn<-read.table(paste0("./",tissue,"/",tissue,"_ensg_hgcn.txt"),header=T,row.names=NULL)
  ensg_hgcn$Name<-gsub("\\..*$", "", ensg_hgcn$Name)
  
  tissue.names<-read_trios(paste0("./",tissue,"/human_",tissue,"_deseq.txt"),inputs=trios_data[tissue])
  nrow(resultsAtoB)
  #bootresults<-cbind(bootresults,tissue.names)
  bootresults$id<-seq(1,nrow(bootresults))
  bootresults$padjust_nomap_not3<-p.adjust(bootresults$p_final,method="BH")
  bootresults$padjust_nomap_t30.1<-p.adjust(ifelse(bootresults$t3<=0.1,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.1<-p.adjust(ifelse(bootresults$t3<=0.1 & bootresults$ismappable==TRUE,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_t30.2<-p.adjust(ifelse(bootresults$t3<=0.2,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.2<-p.adjust(ifelse(bootresults$t3<=0.2 & bootresults$ismappable==TRUE,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_not3wcount<-p.adjust(ifelse(bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_t30.1wcount<-p.adjust(ifelse(bootresults$t3<=0.1 & bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.1wcount<-p.adjust(ifelse(bootresults$t3<=0.1 & bootresults$ismappable==TRUE & bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_t30.2wcount<-p.adjust(ifelse(bootresults$t3<=0.2 & bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.2wcount<-p.adjust(ifelse(bootresults$t3<=0.2 & bootresults$counts>=10 & bootresults$ismappable==TRUE,bootresults$p_final,NA),method="BH")
  bootresults$reverse.pfinal<-resultsBtoA$p_final
  bootresults$reverse.t1<-resultsBtoA$t1
  bootresults$reverse.t2<-resultsBtoA$t2
  bootresults$reverse.t3<-resultsBtoA$t3
  bootresults$reverse.t4<-resultsBtoA$t4
  ensg_hgcn<-read.table(paste0("./",tissue,"/",tissue,"_ensg_hgcn.txt"),header=T,row.names=NULL)
  ensg_hgcn$Name<-gsub("\\..*$", "", ensg_hgcn$Name)
  ensg_hgcn<-ensg_hgcn[!duplicated(ensg_hgcn$Name), ]
  bootresults<-merge(bootresults,ensg_hgcn,by.x='A',by.y='Name',all.x=T)
  colnames(bootresults)[colnames(bootresults) == "Description"] <- "A_hgcn"
  bootresults<-merge(bootresults,ensg_hgcn,by.x='B',by.y='Name',all.x=T)
  colnames(bootresults)[colnames(bootresults) == "Description"] <- "B_hgcn"
  return(bootresults)
}

finaltable_nacond<- function(filterdf,resultsAtoB,resultsBtoA,tissue)
{
  ##union of both ###
  bootresults<-resultsAtoB
  stopifnot(nrow(filterdf)==nrow(resultsAtoB))
  stopifnot(nrow(filterdf)==nrow(resultsBtoA))
  bootresults<-cbind(bootresults,filterdf)
  bootresults$ismappable<-(is.na(bootresults$mapA) | bootresults$mapA>=0.75) & (is.na(bootresults$mapB) | bootresults$mapB>=0.75) & (is.na(bootresults$cmap) | bootresults$cmap==0) 
  bootresults$ismappablewcount<-(is.na(bootresults$mapA) | bootresults$mapA>=0.75) & (is.na(bootresults$mapB) | bootresults$mapB>=0.75) & (is.na(bootresults$cmap) | bootresults$cmap==0) & (bootresults$counts>=10) 
  #converting the gene name form ensg format to hgcn format
  ensg_hgcn<-read.table(paste0("./",tissue,"/",tissue,"_ensg_hgcn.txt"),header=T,row.names=NULL)
  ensg_hgcn$Name<-gsub("\\..*$", "", ensg_hgcn$Name)
  
  tissue.names<-read_trios(paste0("./",tissue,"/human_",tissue,"_deseq.txt"),inputs=trios_data[tissue])
  nrow(resultsAtoB)
  #bootresults<-cbind(bootresults,tissue.names)
  bootresults$id<-seq(1,nrow(bootresults))
  bootresults$padjust_nomap_not3<-p.adjust(bootresults$p_final,method="BH")
  bootresults$padjust_nomap_t30.1<-p.adjust(ifelse(bootresults$t3<=0.1,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.1<-p.adjust(ifelse(bootresults$t3<=0.1 & bootresults$ismappable==TRUE,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_t30.2<-p.adjust(ifelse(bootresults$t3<=0.2,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.2<-p.adjust(ifelse(bootresults$t3<=0.2 & bootresults$ismappable==TRUE,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_not3wcount<-p.adjust(ifelse(bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_t30.1wcount<-p.adjust(ifelse(bootresults$t3<=0.1 & bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.1wcount<-p.adjust(ifelse(bootresults$t3<=0.1 & bootresults$ismappable==TRUE & bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_nomap_t30.2wcount<-p.adjust(ifelse(bootresults$t3<=0.2 & bootresults$counts>=10,bootresults$p_final,NA),method="BH")
  bootresults$padjust_map_t30.2wcount<-p.adjust(ifelse(bootresults$t3<=0.2 & bootresults$counts>=10 & bootresults$ismappable==TRUE,bootresults$p_final,NA),method="BH")
  bootresults$reverse.pfinal<-resultsBtoA$p_final
  bootresults$reverse.t1<-resultsBtoA$t1
  bootresults$reverse.t2<-resultsBtoA$t2
  bootresults$reverse.t3<-resultsBtoA$t3
  bootresults$reverse.t4<-resultsBtoA$t4
  ensg_hgcn<-read.table(paste0("./",tissue,"/",tissue,"_ensg_hgcn.txt"),header=T,row.names=NULL)
  ensg_hgcn$Name<-gsub("\\..*$", "", ensg_hgcn$Name)
  ensg_hgcn<-ensg_hgcn[!duplicated(ensg_hgcn$Name), ]
  bootresults<-merge(bootresults,ensg_hgcn,by.x='A',by.y='Name',all.x=T)
  colnames(bootresults)[colnames(bootresults) == "Description"] <- "A_hgcn"
  bootresults<-merge(bootresults,ensg_hgcn,by.x='B',by.y='Name',all.x=T)
  colnames(bootresults)[colnames(bootresults) == "Description"] <- "B_hgcn"
  return(bootresults)
}


returnmappableindices<-function(summary_gene_map,tissue,summary_cross_map_final)
{
  ## reading the trio names ###
  trionames_map<-read_trios(paste0("./",tissue,"/human_",tissue,"_deseq.txt"),inputs=trios_data[tissue])
  #remove the version of the genes
  trionames_map$A<-sub("\\..*$","",trionames_map$A)
  trionames_map$B<-sub("\\..*$","",trionames_map$B)
  trionames_map$id=1:nrow(trionames_map)
  
  #add gene mappability of A gene 
  trionames_Amerge<-merge(trionames_map,summary_gene_map,by.x='A',by.y='V1')
  trionames_Amerge<-trionames_Amerge[order(trionames_Amerge$id),]
  colnames(trionames_Amerge)<-c("A","L","B","id","A.count","A.max","A.min")
  #add gene mappability of B gene
  trionames_bothmerge<-merge(trionames_Amerge,summary_gene_map,by.x='B',by.y='V1')
  trionames_bothmerge<-trionames_bothmerge[order(trionames_bothmerge$id),]
  colnames(trionames_bothmerge)<-c("B","A","L","id","A.count","A.max","A.min","B.count","B.max","B.min")
  #for gene mappability take the minimum
  trionames_dropped<-trionames_bothmerge[,c("id","L","A","B","A.min","B.min")]
  ## merge the cross mappability 
  trionames_cross_merge<-merge(trionames_dropped,summary_cross_map_final,by.x=c("A","B"),by.y=c("V1","V2"),all.x=T)
  trionames_cross_merge<-trionames_cross_merge[order(trionames_cross_merge$id),]
  #if either of the first gene is 1 then cross map is 0
  trionames_cross_merge_mod<-trionames_cross_merge
  trionames_cross_merge_mod$max_value[trionames_cross_merge_mod$A.min==1 | trionames_cross_merge_mod$B.min==1]<- 0
  trionames_cm_final<-trionames_cross_merge_mod[!is.na(trionames_cross_merge_mod$max_value),]
  #taking the genes with mappability greater than 0.75 and cross mappability 0 
  #trionames_cm_final_filtered<-trionames_cm_final[trionames_cm_final$A.min>=0.75 & trionames_cm_final$B.min>=0.75 & trionames_cm_final$max_value==0 ,]
  
  outdf<-data.frame(mapA=NA,mapB=NA,cmap=NA,counts=NA,stringsAsFactors = FALSE,A=trionames_map$A,B=trionames_map$B)
  outdf$mapA=trionames_cm_final[match(outdf$A,trionames_cm_final$A),"A.min"]
  outdf$mapB=trionames_cm_final[match(outdf$B,trionames_cm_final$B),"B.min"]
  outdf$cmap=trionames_cm_final[match(paste(outdf$A,outdf$B),paste(trionames_cm_final$A,trionames_cm_final$B)),"max_value"]

  ## going to to remove L where the count of unique value is less than 10 ###
  data<-read_data(paste0("./",tissue,"/human_",tissue,"_deseq.txt"),inputs=trios_data[tissue])
  counts=sapply(data,function(x) min(table(x[[1]])))
  outdf$counts<-counts
  return(outdf)
}

prepare_map<-function(x)
{
  #reading the cross mappability file 
  cross_map<-read.table("./mappability/hg38_cross_mappability_strength_symmetric_mean_sorted.txt.gz",sep='\t',header=F)
  #remove the version numbers 
  cross_map$V1<-sub("\\..*$","",cross_map$V1)
  cross_map$V2<-sub("\\..*$","",cross_map$V2)
  #check if any pairs are duplicated 36757 duplicated 
  sum(duplicated(cross_map[,c("V1","V2")]))
  #there are duplicates after removing the versions 
  summary_cross_map<-cross_map %>% group_by(V1,V2) %>% summarize(count=n(),max_value = max(V3),min_value=min(V3))
  
  #create a duplicate and bind it so that we have cis-trans and trans-cis values
  summary_cross_map_rev<-summary_cross_map
  summary_cross_map_rev[,c("V1","V2")]<-summary_cross_map_rev[,c("V2","V1")]
  summary_cross_map_final<-rbind(summary_cross_map,summary_cross_map_rev)
  
  #read the gene mappability file 
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
  save(file="summary_map.rdat",summary_gene_map,summary_cross_map_final)
  
}
# 
# load("summary_map.rdat")
# load("results.muscle.epept.rdat")
# filterdf<-returnmappableindices(summary_gene_map,"muscle",summary_cross_map_final)
# #write.csv(filterdf, file = "filterdf.csv", row.names = FALSE)
# write.csv(filterdf, file = "filterdf_na.csv", row.names = FALSE)
# #bootresults<-finaltable(filterdf,resultsAtoB,resultsBtoA,"muscle")
# bootresults<-finaltable_nacond(filterdf,resultsAtoB,resultsBtoA,"muscle")
# #write.csv(bootresults, file = "bootresults_normal.csv", row.names = FALSE)
# write.csv(bootresults, file = "bootresults_normal_na.csv", row.names = FALSE)
# load("results.muscle.epept.rdat")
# bootresults<-finaltable(filterdf,resultsAtoB,resultsBtoA,"muscle")
# #write.csv(bootresults, file = "bootresults_epept.csv", row.names = FALSE)
# write.csv(bootresults, file = "bootresults_epept_na.csv", row.names = FALSE)


######### prepare tables for each of the tissues ############## 
load("summary_map.rdat")
#pval_type="epept"
prepare_largetable<-function(tissue,pval_type)
{
load(paste0("./datfiles/results.",tissue,".",pval_type,".rdat"))
#tissue<-"muscle"
if(pval_type=="epept")
{
  resultsAtoB <- resultsAtoB[, -which(names(resultsAtoB) %in% c("t3_0", "t3_1","t3_2"))]
  resultsBtoA <- resultsBtoA[, -which(names(resultsBtoA) %in% c("t3_0", "t3_1","t3_2"))]
}
filterdf<-returnmappableindices(summary_gene_map,tissue,summary_cross_map_final)
colnames(resultsAtoB)<-paste0(colnames(resultsAtoB),"_fwd")
colnames(resultsBtoA)<-paste0(colnames(resultsBtoA),"_rev")
temp1<-cbind(filterdf,resultsAtoB,resultsBtoA)

ensg_hgcn<-read.table(paste0("./",tissue,"/",tissue,"_ensg_hgcn.txt"),header=T,row.names=NULL)
ensg_hgcn$Name<-gsub("\\..*$", "", ensg_hgcn$Name)
ensg_hgcn<-ensg_hgcn[!duplicated(ensg_hgcn$Name), ]
temp_Amerge<-merge(temp1,ensg_hgcn,by.x='A',by.y='Name',all.x=T)
colnames(temp_Amerge)[colnames(temp_Amerge) == "Description"] <- "A_hgcn"
temp_Bmerge<-merge(temp_Amerge,ensg_hgcn,by.x='B',by.y='Name',all.x=T)
colnames(temp_Bmerge)[colnames(temp_Bmerge) == "Description"] <- "B_hgcn"

## read gencode v26 gtf file ###
#below 3 lines done
#gencode<-read.table("gencode.v26.GRCh38.genes.genetypeOnly.gtf",sep = ';',stringsAsFactors = FALSE)
#temp<-data.frame(gene=unlist(lapply(gencode$V2, function(x) strsplit(x, " ")[[1]][3])),type=unlist(lapply(gencode$V1, function(x) strsplit(x, " ")[[1]][3])))
#write.table(temp,"gencode.v26.GRCh38.genes.genetypeOnly_processed.gtf",row.names = FALSE)

gencode<-read.table("gencode.v26.GRCh38.genes.genetypeOnly_processed.gtf",header=T,stringsAsFactors = F)
gencode<-gencode[!duplicated(gencode$gene), ]
temp_Amergegen<-merge(temp_Bmerge,gencode,by.x='A_hgcn',by.y='gene',all.x=T)
colnames(temp_Amergegen)[colnames(temp_Amergegen) == "type"] <- "A_type"
temp_Bmergegen<-merge(temp_Amergegen,gencode,by.x='B_hgcn',by.y='gene',all.x=T)
colnames(temp_Bmergegen)[colnames(temp_Bmergegen) == "type"] <- "B_type"
temp_Bmergegen$tissue<-tissue
temp_final<-temp_Bmergegen[,c("tissue","A","B","A_hgcn","B_hgcn","A_type","B_type","mapA","mapB","cmap","counts","p_final_fwd","t1_fwd","t2_fwd","t3_fwd","t4_fwd","oscore_fwd","cseed_fwd","pseed_fwd","p_final_rev","t1_rev","t2_rev","t3_rev","t4_rev","oscore_rev","cseed_rev","pseed_rev")]

return(temp_final);
}
# testing
# normal_table<-data.frame()
# muscle_normal<-prepare_largetable("muscle","normal")
# muslce_epept<-prepare_largetable("muscle","epept")
# adipose_normal<-prepare_largetable("adipose","normal")
# adipose_epept<-prepare_largetable("adipose","epept")
# 
# normal_table<-rbind(normal_table,muscle_normal)
# normal_table<-rbind(normal_table,adipose_normal)
# 
# normal_table<-rbind(normal_table,muscle_normal)
# normal_table<-rbind(normal_table,muslce_epept)


tissues<-c("muscle","Whole_Blood","adipose","Lung","Thyroid","Esophagus_Mucosa","Artery_Tibial","Nerve_Tibial","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg")
normal_table<-data.frame()
for(tissue in tissues)
{
  normal_table<-rbind(normal_table,prepare_largetable(tissue,"normal"))
}
epept_table<-data.frame()
for(tissue in tissues)
{
  epept_table<-rbind(epept_table,prepare_largetable(tissue,"epept"))
}

save(file="final_tables.rdat",normal_table,epept_table)
load("final_tables.rdat")
colnames(normal_table)<-c("Tissue","Gene_A_Source","Gene_B_Target","A_gene_symbol","B_gene_symbol","A_gene_type",
                          "B_gene_type","A_mappability","B_mappability","AB_cross_mappability","min_genotype_group_size",
                          "p_final_forward","Test_1_forward","Test_2_forward","Test_3_forward","Test_4_forward",
                          "Overlap_Score_forward","Child_seed_forward","Parent_seed_forward",
                          "p_final_reverse","Test_1_reverse","Test_2_reverse","Test_3_reverse","Test_4_reverse",
                          "Overlap_Score_reverse","Child_seed_reverse","Parent_seed_reverse")
colnames(epept_table)<-c("Tissue","Gene_A_Source","Gene_B_Target","A_gene_symbol","B_gene_symbol","A_gene_type",
                          "B_gene_type","A_mappability","B_mappability","AB_cross_mappability","min_genotype_group_size",
                          "p_final_forward","Test_1_forward","Test_2_forward","Test_3_forward","Test_4_forward",
                          "Overlap_Score_forward","Child_seed_forward","Parent_seed_forward",
                          "p_final_reverse","Test_1_reverse","Test_2_reverse","Test_3_reverse","Test_4_reverse",
                          "Overlap_Score_reverse","Child_seed_reverse","Parent_seed_reverse")
#a<-muscle_normal
# a<-muslce_epept
# a$filter<- a$mapA>=0.75 & a$mapB>=0.75 & a$cmap==0 & a$counts>=10 & a$t3_fwd<=0.1# & (!grepl("pseudogene",a$A_type)) & (!grepl("pseudogene",a$B_type))
# index<- which(a$filter)
# a$padjust[index]<- p.adjust(a$p_final_fwd[index],method = "BH") 
# a[which(a$padjust<=0.1),]
# 
# write.table(normal_table,"normal_table.csv",sep=',',row.names=F)
# write.table(epept_table,"epept_table.csv",sep=',',row.names=F)
# 
# 
# normal_table_tissue<-normal_table[normal_table$Tissue=='muscle',]
# normal_table_tissue$ismappable<- normal_table_tissue$A_mappability>=0.75 &
#    normal_table_tissue$B_mappability>=0.75 &
#   normal_table_tissue$AB_cross_mappability==0 &
#   normal_table_tissue$min_genotype_group_size>=10 &
#   normal_table_tissue$Test_3_forward<=0.1 
# #normal_table_tissue$ismappable<-(is.na(normal_table_tissue$A_mappability) | normal_table_tissue$A_mappability>=0.75) &
# #                                  (is.na(normal_table_tissue$B_mappability) | normal_table_tissue$B_mappability>=0.75) &
# #                               (is.na(normal_table_tissue$AB_cross_mappability) | normal_table_tissue$AB_cross_mappability==0) &
# #                                   (normal_table_tissue$min_genotype_group_size>=10) &
# #                                normal_table_tissue$Test_3_forward<=0.1 
#                                 
# normal_table_tissue$padj<-p.adjust(ifelse(normal_table_tissue$ismappable==TRUE,normal_table_tissue$p_final_forward,NA),method="BH")
# normal_table_tissue$causal<- ifelse(normal_table_tissue$padj<=0.1 & normal_table_tissue$p_final_reverse>0.1,"causal","NA")
check_for_causal<-function(tissue,big_table,include_na=T)
{
  big_table_tissue<-big_table[big_table$Tissue==tissue,]
  if(include_na==F)
  {
  big_table_tissue$ismappable<- big_table_tissue$A_mappability>=0.75 &
    big_table_tissue$B_mappability>=0.75 &
    big_table_tissue$AB_cross_mappability==0 &
    big_table_tissue$min_genotype_group_size>=10 &
    big_table_tissue$Test_3_forward<=0.1 
  }
  else
  {
  big_table_tissue$ismappable<-(is.na(big_table_tissue$A_mappability) | big_table_tissue$A_mappability>=0.75) &
                                    (is.na(big_table_tissue$B_mappability) | big_table_tissue$B_mappability>=0.75) &
                                 (is.na(big_table_tissue$AB_cross_mappability) | big_table_tissue$AB_cross_mappability==0) &
                                     (big_table_tissue$min_genotype_group_size>=10) &
                                  big_table_tissue$Test_3_forward<=0.1 
  }
  big_table_tissue$padj<-p.adjust(ifelse(big_table_tissue$ismappable==TRUE,big_table_tissue$p_final_forward,NA),method="BH")
  big_table_tissue$causal<- ifelse(big_table_tissue$padj<=0.1 & big_table_tissue$p_final_reverse>0.1,"causal","NA")
  return(big_table_tissue)
}
#checking for muscle
#big_table_muscle<-check_for_causal('muscle',epept_table)
#all tissues
causal_alltissue_na<-data.frame()
causal_alltissue<-data.frame()
alltissues<-c("muscle","Whole_Blood","Skin_Sun_Exposed_Lower_leg","Artery_Tibial","adipose","Thyroid","Nerve_Tibial","Skin_Not_Sun_Exposed_Suprapubic","Lung","Esophagus_Mucosa" )
#doing it for only muscle to get the scatter plot with the updated mappable condition
#alltissues<-c("muscle")
for(tissue in alltissues)
{
big_table_tissue_na<-check_for_causal(tissue,epept_table)
big_table_tissue<-check_for_causal(tissue,epept_table,F)
causal_alltissue_na<-rbind(causal_alltissue_na,big_table_tissue_na[which(big_table_tissue_na$causal=="causal"),])
causal_alltissue<-rbind(causal_alltissue,big_table_tissue[which(big_table_tissue$causal=="causal"),])
}
#View(big_table_tissue)
#save(file="big_table_muscle.rdat",big_table_tissue)
write.table(causal_alltissue_na,"causal_alltissue_na.csv",sep = ",",row.names=F)
write.table(causal_alltissue,"causal_alltissue.csv",sep = ",",row.names=F)

colnames(causal_alltissue)
causal_alltissue_latex<-causal_alltissue[,c("Tissue","Gene_A_Source", "Gene_B_Target" , "A_gene_symbol" ,"B_gene_symbol",
                                            "A_gene_type" ,"B_gene_type","A_mappability" ,"B_mappability" ,"AB_cross_mappability",
                                            "p_final_forward","p_final_reverse"   )]
causal_alltissue_na_latex<-causal_alltissue_na[,c("Tissue","Gene_A_Source", "Gene_B_Target" , "A_gene_symbol" ,"B_gene_symbol",
                                               "A_gene_type" ,"B_gene_type","A_mappability" ,"B_mappability" ,"AB_cross_mappability",
                                               "p_final_forward","p_final_reverse"   )]
write.table(causal_alltissue_na_latex,"causal_alltissue_na_latex.csv",sep = ",",row.names=F)
write.table(causal_alltissue_latex,"causal_alltissue_latex.csv",sep = ",",row.names=F)


