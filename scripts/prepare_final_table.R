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

load("summary_map.rdat")
load("results.muscle.normal.rdat")
filterdf<-returnmappableindices(summary_gene_map,"muscle",summary_cross_map_final)
write.csv(filterdf, file = "filterdf.csv", row.names = FALSE)
bootresults<-finaltable(filterdf,resultsAtoB,resultsBtoA,"muscle")
write.csv(bootresults, file = "bootresults_normal.csv", row.names = FALSE)
load("results.muscle.epept.rdat")
bootresults<-finaltable(filterdf,resultsAtoB,resultsBtoA,"muscle")
write.csv(bootresults, file = "bootresults_epept.csv", row.names = FALSE)





