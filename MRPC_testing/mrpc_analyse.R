
ds<-read.csv("./1000_datasets/MRPC_truth_0_numberofdata_1000_.csv",check.names=FALSE,header=TRUE)
read.csv("./1000_datasets/MRPC_truth_1_numberofdata_1000_.csv")
read.csv("./1000_datasets/PC_truth_0_numberofdata_1000_.csv")
read.csv("./1000_datasets/PC_truth_1_numberofdata_1000_.csv")


#use unique to find the number of configurations which turned up in the result
#unique(ds)

temp=vector("character",1000)
ds$model<-apply(ds,1,function(x) paste0(unlist(labels[which(x %in% 1)]),collapse=","))
