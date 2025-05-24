library(ggplot2)
dist_point_line <- function(a, slope, intercept) {
  b = c(1, intercept+slope)
  c = c(-intercept/slope,0)       
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  return(abs(det(m))/sqrt(sum(v1*v1)))
}

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



### causal data ######## 


causal_data<- read.table("./mpmicorspear_causal_wilko1752.csv", header = TRUE, sep = ",", dec = ".")
causal_data$abs_spear_causal<- abs(causal_data$spear_causal)

obj<- lm(causal_data$abs_spear_causal ~ causal_data$bcmi_causal)
c<-coef(obj)[1]
m<-coef(obj)[2]
summary(obj)

fitted <- obj$fitted.values
smaller <-  causal_data$abs_spear_causal <= fitted
causal_data$lowerside<- causal_data$abs_spear_causal <= fitted
causal_data$distance<-apply(causal_data,1,function(x) dist_point_line(c(x['bcmi_causal'],x['abs_spear_causal']), slope = m, intercept = c) )


### independent data ###### 

indep_data<- read.table("./mpmicorspear_indep_wilko1752.csv", header = TRUE, sep = ",", dec = ".")
indep_data$abs_spear_indep<- abs(indep_data$spear_indep)
indep_data$diff<- abs(indep_data$bcmi_indep - indep_data$abs_spear_indep)
obj<- lm(indep_data$abs_spear_indep ~ indep_data$bcmi_indep)
c<-coef(obj)[1]
m<-coef(obj)[2]
summary(obj)

fitted <- obj$fitted.values
smaller <-indep_data$abs_spear_indep <= fitted
indep_data$lowerside<- indep_data$abs_spear_indep <= fitted

indep_data$distance<-apply(indep_data,1,function(x) dist_point_line(c(x['bcmi_indep'],x['abs_spear_indep']), slope = m, intercept = c) )


##########################



##### going to try 10 different seeds but with cut off point different for spearman############
### the cutoff point is corresponding to the bcmi cut off ### 
######### sir suggestion ##############
base.dir="./mpmi_indices_1752/"
seed=sample(.Machine$integer.max, 1) # 1140350788
set.seed(1140350788)
num_seeds<-10
random_seeds <- sample.int(2^31 - 1, num_seeds)
# 1588846348 1463234807 1356651008  704118021 1696680727  671479058 1703830755  900040367  175283650  1732477378
### bcmi > 0.05 0.1 0.15 ##########
seed_count=0
for (s in random_seeds)
{
  seed_count=seed_count+1
  base.dir=paste0("./mpmi_indices_1752/seed",seed_count,"/")
  dir.create(base.dir)
  print(s)
  set.seed(s) 
  ## bcmi
  for(l in c(0.05,0.1,0.15))
  {
    causal_points=c()
    for(i in c(0,0.05,0.1))
    {
      points=which(  causal_data$lowerside==TRUE & causal_data$distance>i & causal_data$bcmi_causal > l) #new condition added
      points=points-1 # to make them zero indexed for python 
      lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_bcmi",l,"_",i,".csv"),append=TRUE)
      causal_points=c(causal_points,length(points))
      cat("bcmi",l,"odc cutoff",i," number of points", length(points),"\n")
    }
    cutoff=c(0,0.05,0.1)
    j=0
    for(i in causal_points)
    {
      j=j+1
      if( i > sum(indep_data$bcmi_indep>l) )
      {
        print("independent data doesnt have enought points so taking all of them")
        #print(sum(indep_data$bcmi_indep>l))
        #print(l)
        #print(i)
        if(sum(indep_data$bcmi_indep>l)==0)
        {
          print("number of points in independent data set is zero")
        }
        points=which(indep_data$bcmi_indep>l)
      }
      else
      {
        points=sample(which(indep_data$bcmi_indep>l),size=i,replace = FALSE) #new condition added
      }
      
      points=points-1 # to make them zero indexed for python 
      lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_bcmi",l,"_",cutoff[j],".csv"),append=TRUE)
    }
    
  }
  
  #### spearman 0.05 0.1 0.15 ### 
  
  for(l in c(0.05,0.1,0.15))
  {
    causal_points=c()
    obj<- lm(abs_spear_causal ~ bcmi_causal,data=causal_data)
    v=predict(obj,data.frame(bcmi_causal=c(l)))
    for(i in c(0,0.05,0.1))
    {
      points=which(   causal_data$lowerside==TRUE & causal_data$distance>i & causal_data$abs_spear_causal > v) #new condition added
      points=points-1 # to make them zero indexed for python 
      lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_spear",l,"_",i,".csv"),append=TRUE)
      causal_points=c(causal_points,length(points))
      cat("spear",l,"odc cutoff",i," number of points", length(points),"\n")
      
    }
    cutoff=c(0,0.05,0.1)
    j=0
    #obj<- lm(abs_spear_indep ~ bcmi_indep,data=indep_data)
    #v=predict(obj,data.frame(bcmi_indep=c(l)))
    for(i in causal_points)
    {
      j=j+1
      if( i > sum(indep_data$abs_spear_indep>v) )
      {
        print("independent data doesnt have enought points so taking all of them")
        if(sum(indep_data$abs_spear_indep>l)==0)
        {
          print("number of points in independent data set is zero")
        }
        points=which(indep_data$abs_spear_indep>v)
      }
      else
      {
        points=sample(which(indep_data$abs_spear_indep>v),size=i,replace = FALSE) #new condition added
      }
      
      points=points-1 # to make them zero indexed for python 
      lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_spear",l,"_",cutoff[j],".csv"),append=TRUE)
    }
    
  }
  
  ### spearman or bcmi 0.05 0.1 0.15 ###### 
  spearman_cutoff_causal=c()
  spearman_cutoff_indep=c()
  for(l in c(0.05,0.1,0.15))
  {
    
    causal_points=c()
    obj<- lm(abs_spear_causal ~ bcmi_causal,data=causal_data)
    v=predict(obj,data.frame(bcmi_causal=c(l)))
    spearman_cutoff_causal=c(spearman_cutoff_causal,v)
    for(i in c(0,0.05,0.1))
    {
      
      points=which(causal_data$lowerside==TRUE & causal_data$distance>i & (causal_data$bcmi_causal > l  | causal_data$abs_spear_causal>v)) #new condition added
      points=points-1 # to make them zero indexed for python 
      lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_both",l,"_",i,".csv"),append=TRUE)
      causal_points=c(causal_points,length(points))
      cat("both bcmi",l,"spearman",v, "odc cutoff",i," number of points", length(points),"\n")
    }
    cutoff=c(0,0.05,0.1)
    j=0
    #obj<- lm(abs_spear_indep ~ bcmi_indep,data=indep_data)
    #v=predict(obj,data.frame(bcmi_indep=c(l)))
    #spearman_cutoff_indep=c(spearman_cutoff_indep,v)
    for(i in causal_points)
    {
      j=j+1
      if( i > sum(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>v) )
      {
        print("independent data doesnt have enought points so taking all of them")
        if(sum(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>v) ==0)
        {
          print("number of points in independent data set is zero")
        }
        points=which(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>v)
      }
      else
      {
        points=sample(which(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>v),size=i,replace = FALSE) #new condition added
      }
      
      points=points-1 # to make them zero indexed for python 
      lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_both",l,"_",cutoff[j],".csv"),append=TRUE)
    }
    
  }
  
}

write.table(spearman_cutoff_causal,file="./mpmi_indices_1752/spearman_cutoff_causal.txt",quote=FALSE,row.names=F,col.names=F)
write.table(spearman_cutoff_indep,file="./mpmi_indices_1752/spearman_cutoff_indep.txt",quote=FALSE,row.names=F,col.names=F)


########### Figure 4a and S8 ##############
library(gridExtra)
library(grid)
### causal data Figure 4 a
nlcdyeast_causal<-read.table("./yeast_nlcd_result/yeast_causal_1752.csv",header=T,sep=',') 
nlcdyeast_indp<-read.table("./yeast_nlcd_result/yeast_indp_1752.csv",header=T,sep=',') 
cityeast_causal<-read.table("./yeast_cit_result/yeast_causal_1752.csv",header=T,sep=',') 
cityeast_indp<-read.table("./yeast_cit_result/yeast_indp_1752.csv",header=T,sep=',')
nlcd_plot<-ggplot(causal_data, aes(x=bcmi_causal, y=abs_spear_causal, color=as.factor(ifelse(nlcdyeast_causal$p_final < 0.05, "Causal (p<=0.05)", "Independent (p>0.05)")))) + 
  geom_point(alpha=0.5)+
  geom_smooth(method='lm', color = 'black')+ggtitle("NLCD ",) + labs(x="BCMI",y=expression(paste("Absolute Spearman (|", rho, "|)", sep = "")),colour = "NLCD prediction")+theme_bw()+ theme(legend.position="none",plot.title = element_text(size=10,hjust = 0.5)) 
cit_plot<-ggplot(causal_data, aes(x=bcmi_causal, y=abs_spear_causal, color=as.factor(ifelse(cityeast_causal$p_cit < 0.05, "Causal (p<=0.05)", "Independent (p>0.05)")))) + 
  geom_point(alpha=0.5)+
  geom_smooth(method='lm', color = 'black')+ggtitle("CIT ") +labs(x="BCMI",colour = "NLCD/CIT predictions:") +ylab(NULL)+ theme_bw()+ theme(legend.position = "bottom",axis.text.y=element_blank(),plot.title = element_text(size=10,hjust = 0.5))
#legend extractor function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(cit_plot)

combined_plot <- grid.arrange(arrangeGrob(nlcd_plot + theme(legend.position="none"),
                                          cit_plot + theme(legend.position="none"),
                                          nrow=1,widths = c(1.1, 0.96)),
                              mylegend, nrow=2,heights=c(10, 1))

## going to do the same for independent data ##  Figure S8
nlcd_plot<-ggplot(indep_data, aes(x=bcmi_indep, y=abs_spear_indep, color=as.factor(ifelse(nlcdyeast_indp$p_final < 0.05, "Causal (p<=0.05)", "Independent (p>0.05)")))) + 
  geom_point(alpha=0.5)+
  geom_smooth(method='lm', color = 'black')+ggtitle("NLCD ",) + labs(x="BCMI",y=expression(paste("Absolute Spearman (|", rho, "|)", sep = "")),colour = "NLCD prediction")+theme_bw()+ theme(legend.position="none",plot.title = element_text(size=10,hjust = 0.5)) 
cit_plot<-ggplot(indep_data, aes(x=bcmi_indep, y=abs_spear_indep, color=as.factor(ifelse(cityeast_indp$p_cit < 0.05, "Causal (p<=0.05)", "Independent (p>0.05)")))) + 
  geom_point(alpha=0.5)+
  geom_smooth(method='lm', color = 'black')+ggtitle("CIT ") +labs(x="BCMI",colour = "NLCD/CIT predictions:") +ylab(NULL)+ theme_bw()+ theme(legend.position = "bottom",axis.text.y=element_blank(),plot.title = element_text(size=10,hjust = 0.5))
#legend extractor function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(cit_plot)

combined_plot <- grid.arrange(arrangeGrob(nlcd_plot + theme(legend.position="none"),
                                          cit_plot + theme(legend.position="none"),
                                          nrow=1,widths = c(1.1, 0.94)),
                              mylegend, nrow=2,heights=c(10, 1)) #c(10,1)
######################################

#saving files for the supplementary file 
supp_data_causal<- causal_data[,c('bcmi_causal','abs_spear_causal','distance')]
supp_data_indep<-indep_data[,c('bcmi_indep','abs_spear_indep','distance')]
write.table(supp_data_causal,"supp_causal_yeast_1752.csv",quote = FALSE,row.names = FALSE)
write.table(supp_data_indep,"supp_indep_yeast_1752.csv",quote = FALSE,row.names = FALSE)
