tissues<-c("muscle","Whole_Blood","adipose","Lung","Thyroid","Esophagus_Mucosa","Artery_Tibial","Nerve_Tibial","Skin_Not_Sun_Exposed_Suprapubic","Skin_Sun_Exposed_Lower_leg")
for(tissue in tissues)
{
#tissue<-"muscle"
t1pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/pvalues_t1loss_noheader_",tissue_name[tissue],".txt"),sep = ''))
t2pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/pvalues_t2loss_noheader_",tissue_name[tissue],".txt"),sep = ''))
t3_0pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/pvalues_t3loss_0_noheader_",tissue_name[tissue],".txt"),sep = ''))
t3_1pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/pvalues_t3loss_1_noheader_",tissue_name[tissue],".txt"),sep = ''))
t3_2pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/pvalues_t3loss_2_noheader_",tissue_name[tissue],".txt"),sep = ''))
t4pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/pvalues_t4loss_noheader_",tissue_name[tissue],".txt"),sep = ''))

ptable<-data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)
ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,ptable$t3_2,na.rm = T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)
resultsAtoB<-ptable

t1pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/pvalues_t1loss_noheader_",tissue_name[tissue],".txt"),sep = ''))
t2pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/pvalues_t2loss_noheader_",tissue_name[tissue],".txt"),sep = ''))
t3_0pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/pvalues_t3loss_0_noheader_",tissue_name[tissue],".txt"),sep = ''))
t3_1pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/pvalues_t3loss_1_noheader_",tissue_name[tissue],".txt"),sep = ''))
t3_2pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/pvalues_t3loss_2_noheader_",tissue_name[tissue],".txt"),sep = ''))
t4pval<-t(read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/pvalues_t4loss_noheader_",tissue_name[tissue],".txt"),sep = ''))

ptable<-data.frame(t1=t1pval,t2=t2pval,t3_0=t3_0pval,t3_1=t3_1pval,t3_2=t3_2pval,t4=t4pval)
ptable$t3<-pmin(ptable$t3_0,ptable$t3_1,ptable$t3_2,na.rm = T)
ptable$p_final<-pmax(ptable$t1,ptable$t2,ptable$t3,ptable$t4)
resultsBtoA<-ptable

### transfer the seeds and oscore to epept from the normal table ###
resultsAtoB_normal<-read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/normal/test_",tissue_name[tissue],".csv"),sep =',',header=T,numerals = "no.loss",stringsAsFactors = F)
resultsBtoA_normal<-read.table(paste0("./../results/final_human/humans_tissue/",tissue,"/reverse/test_",tissue_name[tissue],"_rev.csv"),sep =',',header=T,numerals = "no.loss",stringsAsFactors = F)
resultsAtoB_normal$parent_seed<-resultsAtoB_normal[1,"parent_seed"]
resultsBtoA_normal$parent_seed<-resultsBtoA_normal[1,"parent_seed"]
colnames(resultsAtoB_normal)<-c("p_final","t1","t2","t3","t4","oscore","cseed","pseed")
colnames(resultsBtoA_normal)<-c("p_final","t1","t2","t3","t4","oscore","cseed","pseed")
resultsAtoB$pseed<-resultsAtoB_normal$pseed
resultsAtoB$cseed<-resultsAtoB_normal$cseed
resultsAtoB$oscore<-resultsAtoB_normal$oscore
resultsBtoA$pseed<-resultsBtoA_normal$pseed
resultsBtoA$cseed<-resultsBtoA_normal$cseed
resultsBtoA$oscore<-resultsBtoA_normal$oscore

save(file=paste0("./datfiles/results.",tissue,".epept.rdat"),resultsAtoB,resultsBtoA) 
#reassign the names
resultsAtoB<-resultsAtoB_normal
resultsBtoA<-resultsBtoA_normal
save(file=paste0("./datfiles/results.",tissue,".normal.rdat"),resultsAtoB,resultsBtoA) 

}





