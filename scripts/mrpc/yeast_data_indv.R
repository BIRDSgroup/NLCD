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

library(MRPC)
types<-c("1","0")
for (type in types)
{
    yeast_data<-read_data(paste0("/home/aravind/Documents/MRPC/yeast_data/yeast_trios/yeastgt_",type,"_wilko1752_ready.txt"),inputs=1752) #100 is the number of trios 
    res_values<-list()
    res_list<-list()
    res_trio<-list()
    counter=0;
    for(i in 1:1752)
    {
      print(paste0("trio number ",i))
      for(f in c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1))
      {
        counter=counter+1;
        res_trio[[counter]]<-c(i,f)
        yeast_df<-data.frame(L=unlist(yeast_data[[i]][1]),A=unlist(yeast_data[[i]][2]),B=unlist(yeast_data[[i]][3]))
        
        n <- nrow(yeast_df)
        # Extract the node/column names
        V <- colnames(yeast_df)
        # Calculate Pearson correlation
        #suffStat_C <- list(C = cor(sim_df), n = n)
        # Calculate Robust correlation 
        robustCor_obj<-RobustCor(yeast_df,Beta=0.005,plot=FALSE)
        suffStat_C<-list(C=robustCor_obj$RR,n=n)
        
        # Infer the graph by MRPC
        yeast_df.mrpc <- MRPC(data = yeast_df,
                            suffStat = suffStat_C,
                            GV = 1,
                            FDR = f,
                            indepTest = 'gaussCItest',
                            labels = V,
                            FDRcontrol = 'LOND',
                            verbose = FALSE)
        
        
        #res_list[[i]]<-sim_df.mrpc@pval
        #res_list[[i]]<-names(sim_df.mrpc@graph@edgeData)
        res_list[[counter]]<-names(yeast_df.mrpc@graph@edgeData@data)
        res_values[[counter]]<-unname(unlist(yeast_df.mrpc@graph@edgeData@data))
      }
    }
    max_length <- max(sapply(res_list, length))
    
    list_of_vectors_names<-lapply(res_list,function(x){
      if(length(x)==0)
        x<-"NA"
      length(x) <- max_length
      x
    })
    res_table_names <- do.call(rbind,list_of_vectors_names)
    
    list_of_vectors_values<-lapply(res_values,function(x){
      if(length(x)==0)
        x<-"NA"
      length(x) <- max_length
      x
    })
    res_table_values <- do.call(rbind,list_of_vectors_values)
    
    res_table_trios<-do.call(rbind,res_trio)
    
    colnames(res_table_values)<-rep("weight values",dim(res_table_values)[2])
    colnames(res_table_names)<-rep("edge names",dim(res_table_names)[2])
    colnames(res_table_trios)<-c("trio number","fdr")
    
    
    res_table<-cbind(res_table_trios,res_table_names,res_table_values)
    
    
    write.csv(res_table,file=paste0("/home/aravind/Documents/MRPC/yeast_results/yeast_indv/",type,"_trios_mrpc",".csv"),row.names = F,quote = F)
    print(paste0(type," completed"))

}

