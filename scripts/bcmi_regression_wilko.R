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
yeast_data<-read_data('yeastgt_1_wilko1234_ready.txt',input=1234)
causal_data<- read.table("./mpmicorspear_causal_wilko1234.csv", header = TRUE, sep = ",", dec = ".")
ggplot(causal_data,aes(mi_causal,bcmi_causal))+geom_point()+geom_abline(color='red')+ggtitle("Scatter plot of bcmi vs spearman causal data")
ggsave("miVSbcmi_causal_wilko1234.png")
causal_data$abs_spear_causal<- abs(causal_data$spear_causal)

obj<- lm(causal_data$abs_spear_causal ~ causal_data$bcmi_causal)
c<-coef(obj)[1]
m<-coef(obj)[2]
summary(obj)
#ggplot(causal_data, aes(x =bcmi_causal, y =abs_spear_causal)) +geom_point() +stat_smooth(method = "lm", se=FALSE)

#function to caluclate the distance from a line given slope and intercept to a point 

fitted <- obj$fitted.values
smaller <-      causal_data$abs_spear_causal <= fitted
causal_data$lowerside<- causal_data$abs_spear_causal <= fitted
# plot of regression line and the points 
ggplot(causal_data, aes(x=bcmi_causal, y=abs_spear_causal, color=as.factor(lowerside))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+
  scale_color_manual(values = c("green", "yellow"))
ggsave("regressionplot_wilko1234_causal.png")

causal_data$distance<-apply(causal_data,1,function(x) dist_point_line(c(x['bcmi_causal'],x['abs_spear_causal']), slope = m, intercept = c) )
#table(causal_data$lowerside==TRUE & causal_data$distance>0.01) old values
table(causal_data$lowerside==TRUE & causal_data$distance>0.1 & causal_data$bcmi_causal > 0.05)
#table(causal_data$distance>0 & causal_data$bcmi_causal > 0.05)
#table( causal_data$distance>0)
point1indices<- which(causal_data$lowerside==TRUE & causal_data$distance>0.15)


seed=sample(.Machine$integer.max, 1) # 1140350788
set.seed(1140350788)
causal_points=c()
for(i in seq(0,0.1,by=0.01))
{
  points=which( causal_data$lowerside==TRUE & causal_data$distance>i & causal_data$bcmi_causal > 0.05) #new condition added
  points=points-1 # to make them zero indexed for python 
  lapply(points, write, paste0("yeast_wilko_mi_causal_",i,".csv"),append=TRUE)
  causal_points=c(causal_points,length(points))
  #points=sample(seq(1,1234,by=1),size=length(points),replace=FALSE)
  #points=points-1
  #lapply(points, write, paste0("yeast_wilko_mi_indep_",i,".csv"),append=TRUE)
  
}
# 28430 31414 48803 > 0.2 
# 255 673 4290 >0.1
for (i in point1indices){
  L=unlist(yeast_data[[i]][1])
  A=unlist(yeast_data[[i]][2])
  B=unlist(yeast_data[[i]][3])
  highdf<-data.frame(L,A,B)
  ggplot(highdf,aes(A,B,color=L))+geom_point()+stat_smooth(method="loess")+ggtitle(paste0("Orthogonal distance >0.15 index= ",i))
  ggsave(paste0('scattertrio_wilko1234',i,'causal0.15.png'))
}

ggplot(causal_data,aes(bcmi_causal,abs_spear_causal,color=distance>0.15))+geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("distance from regress line cut off=0.15")
ggsave('mispear_causal_wilko1234.png')
#calculate the distance from each of these points 

###independent data######

indep_data<- read.table("./mpmicorspear_indep_wilko1234.csv", header = TRUE, sep = ",", dec = ".")
ggplot(indep_data,aes(mi_indep,bcmi_indep))+geom_point()+geom_abline(color='red')+ggtitle("Scatter plot of bcmi vs spearman independent data")
ggsave("miVSbcmi_indep_wilko1234.png")
indep_data$abs_spear_indep<- abs(indep_data$spear_indep)
indep_data$diff<- abs(indep_data$bcmi_indep - indep_data$abs_spear_indep)
obj<- lm(indep_data$abs_spear_indep ~ indep_data$bcmi_indep)
c<-coef(obj)[1]
m<-coef(obj)[2]
summary(obj)
#ggplot(indp_data, aes(x =bcmi_indp, y =abs_spear_indp)) +geom_point() +stat_smooth(method = "lm", se=FALSE)

#function to caluclate the distance from a line given slope and intercept to a point 
dist_point_line <- function(a, slope, intercept) {
  b = c(1, intercept+slope)
  c = c(-intercept/slope,0)       
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  return(abs(det(m))/sqrt(sum(v1*v1)))
}

fitted <- obj$fitted.values
smaller <-indep_data$abs_spear_indep <= fitted
indep_data$lowerside<- indep_data$abs_spear_indep <= fitted
# plot of regression line and the points 
ggplot(indep_data, aes(x=bcmi_indep, y=abs_spear_indep, color=as.factor(smaller))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+
  scale_color_manual(values = c("green", "yellow"))
ggsave("regressionplot_indep_wilko1234.png")

indep_data$distance<-apply(indep_data,1,function(x) dist_point_line(c(x['bcmi_indep'],x['abs_spear_indep']), slope = m, intercept = c) )
#table(indep_data$lowerside==TRUE & indep_data$distance<0.1)
table(indep_data$bcmi_indep>0.05)
ggplot(indep_data,aes(bcmi_indep,abs_spear_indep,color=distance>0.15))+geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("distance from regress line cut off=0.15")
ggsave('mispear_indp_wilko1234.png')
#calculate the distance from each of these points 
yeast_data<-read_data('./../findr/findrfiles/scripts/yeastgt_0_wilko1234_ready.txt',input=1234)
point1indices<- which(indp_data$lowerside==TRUE & indp_data$distance>0.15)
# 28430 31414 48803 > 0.2 
# 255 673 4290 >0.1
for(i in causal_points)
{
  print(i)
  points=sample(which(indep_data$bcmi_indep>0.05),size=i,replace = FALSE) #new condition added
  
  points=points-1 # to make them zero indexed for python 
  lapply(points, write, paste0("yeast_wilko_mi_indep_",i,".csv"),append=TRUE)
  #points=sample(seq(1,1234,by=1),size=length(points),replace=FALSE)
  #points=points-1
  #lapply(points, write, paste0("yeast_wilko_mi_indep_",i,".csv"),append=TRUE)
  
}
for (i in point1indices){
  L=unlist(yeast_data[[i]][1])
  A=unlist(yeast_data[[i]][2])
  B=unlist(yeast_data[[i]][3])
  highdf<-data.frame(L,A,B)
  ggplot(highdf,aes(A,B,color=L))+geom_point()+ggtitle(paste0("Scatter plot of trio where the orthogonal distance >0.15 index= ",i))
  ggsave(paste0('scattertrio_wilko',i,'indep0.15.png'))
}

######## plotting the p values on the graphs ##########
nlcdyeast_causal<-read.table("yeast_causal1.csv",header=T,sep=',') 
nlcdyeast_indp<-read.table("yeast_indp1.csv",header=T,sep=',') 
cityeast_causal<-read.table("yeast_causal_1.csv",header=T,sep=',') 
cityeast_indp<-read.table("yeast_indp_1.csv",header=T,sep=',') 

nlcdyeast_causal$p_final

ggplot(causal_data, aes(x=bcmi_causal, y=abs_spear_causal, color=as.factor(ifelse(nlcdyeast_causal$p_final < 0.05, "Causal (p<=0.05)", "Independent (p>0.05")))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("NLCD p-value yeast causal",)+ labs(x="Bias corrected Mutual Information (BCMI)",y="Absolute Spearman",colour = "NLCD prediction")+ theme_bw() 
#changing png to eps 
ggsave("nlcd_yeast_causal.eps",dpi = 1200)
svg("nlcd_yeast_causal.svg")
ggplot(causal_data, aes(x=bcmi_causal, y=abs_spear_causal, color=as.factor(ifelse(cityeast_causal$p_cit < 0.05, "Causal (p<=0.05)", "Independent (p>0.05")))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("CIT p-value yeast causal")+ labs(x="Bias corrected Mutual Information (BCMI)",y="Absolute Spearman",colour = "CIT prediction") + theme_bw()
ggsave("cit_yeast_causal.eps",dpi = 1200)
ggplot(indep_data, aes(x=bcmi_indep, y=abs_spear_indep, color=as.factor(ifelse(nlcdyeast_indp$p_final < 0.05, "Causal (p<=0.05)", "Independent (p>0.05")))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("NLCD p-value yeast independent") + labs(x="Bias corrected Mutual Information (BCMI)",y="Absolute Spearman",colour = "NLCD prediction")+theme_bw()
ggsave("nlcd_yeast_indp.eps",dpi = 1200)
ggplot(indep_data, aes(x=bcmi_indep, y=abs_spear_indep, color=as.factor(ifelse(cityeast_indp$p_cit < 0.05, "Causal (p<=0.05)", "Independent (p>0.05")))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("CIT p-value yeast independent") + labs(x="Bias corrected Mutual Information (BCMI)",y="Absolute Spearman",colour = "CIT prediction") +theme_bw()
ggsave("cit_yeast_indp.eps",dpi = 1200)
######## analysis for below, top and all ############## 

### causal data ######## 


causal_data<- read.table("./mpmicorspear_causal_wilko1234.csv", header = TRUE, sep = ",", dec = ".")
ggplot(causal_data,aes(mi_causal,bcmi_causal))+geom_point()+geom_abline(color='red')+ggtitle("Scatter plot of bcmi vs spearman causal data")
causal_data$abs_spear_causal<- abs(causal_data$spear_causal)

obj<- lm(causal_data$abs_spear_causal ~ causal_data$bcmi_causal)
c<-coef(obj)[1]
m<-coef(obj)[2]
summary(obj)

fitted <- obj$fitted.values
smaller <-      causal_data$abs_spear_causal <= fitted
causal_data$lowerside<- causal_data$abs_spear_causal <= fitted
# plot of regression line and the points 
ggplot(causal_data, aes(x=bcmi_causal, y=abs_spear_causal, color=as.factor(lowerside))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+
  scale_color_manual(values = c("green", "yellow"))

causal_data$distance<-apply(causal_data,1,function(x) dist_point_line(c(x['bcmi_causal'],x['abs_spear_causal']), slope = m, intercept = c) )
#table(causal_data$lowerside==TRUE & causal_data$distance>0.01) old values
table(causal_data$lowerside==TRUE & causal_data$distance>0.1 & causal_data$bcmi_causal > 0.05)

### independent data ###### 

indep_data<- read.table("./mpmicorspear_indep_wilko1234.csv", header = TRUE, sep = ",", dec = ".")
ggplot(indep_data,aes(mi_indep,bcmi_indep))+geom_point()+geom_abline(color='red')+ggtitle("Scatter plot of bcmi vs spearman independent data")
indep_data$abs_spear_indep<- abs(indep_data$spear_indep)
indep_data$diff<- abs(indep_data$bcmi_indep - indep_data$abs_spear_indep)
obj<- lm(indep_data$abs_spear_indep ~ indep_data$bcmi_indep)
c<-coef(obj)[1]
m<-coef(obj)[2]
summary(obj)

fitted <- obj$fitted.values
smaller <-indep_data$abs_spear_indep <= fitted
indep_data$lowerside<- indep_data$abs_spear_indep <= fitted
# plot of regression line and the points 
ggplot(indep_data, aes(x=bcmi_indep, y=abs_spear_indep, color=as.factor(smaller))) + 
  geom_point()+
  geom_smooth(method='lm', color = 'black')+
  scale_color_manual(values = c("green", "yellow"))

indep_data$distance<-apply(indep_data,1,function(x) dist_point_line(c(x['bcmi_indep'],x['abs_spear_indep']), slope = m, intercept = c) )
#table(indep_data$lowerside==TRUE & indep_data$distance<0.1)
table(indep_data$bcmi_indep>0.05)
ggplot(indep_data,aes(bcmi_indep,abs_spear_indep,color=distance>0.15))+geom_point()+
  geom_smooth(method='lm', color = 'black')+ggtitle("distance from regress line cut off=0.15")


######### sir suggestion ##############
base.dir="./mpmi_indices/"
seed=sample(.Machine$integer.max, 1) # 1140350788
set.seed(1140350788)
### bcmi > 0.05 0.1 0.15 ##########
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
  for(i in c(0,0.05,0.1))
  {
    points=which(   causal_data$lowerside==TRUE & causal_data$distance>i & causal_data$abs_spear_causal > l) #new condition added
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_spear",l,"_",i,".csv"),append=TRUE)
    causal_points=c(causal_points,length(points))
    cat("spear",l,"odc cutoff",i," number of points", length(points),"\n")
    
  }
  cutoff=c(0,0.05,0.1)
  j=0
  for(i in causal_points)
  {
    j=j+1
    if( i > sum(indep_data$abs_spear_indep>l) )
    {
      print("independent data doesnt have enought points so taking all of them")
      points=which(indep_data$abs_spear_indep>l)
    }
    else
    {
      points=sample(which(indep_data$abs_spear_indep>l),size=i,replace = FALSE) #new condition added
    }
    
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_spear",l,"_",cutoff[j],".csv"),append=TRUE)
  }
  
}

### spearman or bcmi 0.05 0.1 0.15 ###### 
for(l in c(0.05,0.1,0.15))
{
  causal_points=c()
  for(i in c(0,0.05,0.1))
  {
    points=which(  causal_data$lowerside==TRUE & causal_data$distance>i & (causal_data$bcmi_causal > l  | causal_data$abs_spear_causal>l)) #new condition added
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_both",l,"_",i,".csv"),append=TRUE)
    causal_points=c(causal_points,length(points))
    cat("both",l,"odc cutoff",i," number of points", length(points),"\n")
  }
  cutoff=c(0,0.05,0.1)
  j=0
  for(i in causal_points)
  {
    j=j+1
    if( i > sum(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>l) )
    {
      print("independent data doesnt have enought points so taking all of them")
      points=which(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>l)
    }
    else
    {
      points=sample(which(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>l),size=i,replace = FALSE) #new condition added
    }
    
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_both",l,"_",cutoff[j],".csv"),append=TRUE)
  }
  
}

##### going to try 10 different seeds ############
######### sir suggestion ##############
base.dir="./mpmi_indices/"
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
  base.dir=paste0("./mpmi_indices/seed",seed_count,"/")
  dir.create(base.dir)
  print(s)
  set.seed(s) 
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
  for(i in c(0,0.05,0.1))
  {
    points=which(   causal_data$lowerside==TRUE & causal_data$distance>i & causal_data$abs_spear_causal > l) #new condition added
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_spear",l,"_",i,".csv"),append=TRUE)
    causal_points=c(causal_points,length(points))
    cat("spear",l,"odc cutoff",i," number of points", length(points),"\n")
    
  }
  cutoff=c(0,0.05,0.1)
  j=0
  for(i in causal_points)
  {
    j=j+1
    if( i > sum(indep_data$abs_spear_indep>l) )
    {
      print("independent data doesnt have enought points so taking all of them")
      points=which(indep_data$abs_spear_indep>l)
    }
    else
    {
      points=sample(which(indep_data$abs_spear_indep>l),size=i,replace = FALSE) #new condition added
    }
    
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_spear",l,"_",cutoff[j],".csv"),append=TRUE)
  }
  
}

### spearman or bcmi 0.05 0.1 0.15 ###### 
for(l in c(0.05,0.1,0.15))
{
  causal_points=c()
  for(i in c(0,0.05,0.1))
  {
    points=which(  causal_data$lowerside==TRUE & causal_data$distance>i & (causal_data$bcmi_causal > l  | causal_data$abs_spear_causal>l)) #new condition added
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_causal_both",l,"_",i,".csv"),append=TRUE)
    causal_points=c(causal_points,length(points))
    cat("both",l,"odc cutoff",i," number of points", length(points),"\n")
  }
  cutoff=c(0,0.05,0.1)
  j=0
  for(i in causal_points)
  {
    j=j+1
    if( i > sum(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>l) )
    {
      print("independent data doesnt have enought points so taking all of them")
      points=which(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>l)
    }
    else
    {
      points=sample(which(indep_data$bcmi_indep>l | indep_data$abs_spear_indep>l),size=i,replace = FALSE) #new condition added
    }
    
    points=points-1 # to make them zero indexed for python 
    lapply(points, write, paste0(base.dir,"yeast_wilko_mi_indep_both",l,"_",cutoff[j],".csv"),append=TRUE)
  }
  
}
  
}


##### going to try 10 different seeds but with cut off point different for spearman############
### the cutoff point is corresponding to the bcmi cut off ### 
######### sir suggestion ##############
base.dir="./mpmi_indices/"
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
  base.dir=paste0("./mpmi_indices/seed",seed_count,"/")
  dir.create(base.dir)
  print(s)
  set.seed(s) 
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

write.table(spearman_cutoff_causal,file="spearman_cutoff_causal.txt",quote=FALSE,row.names=F,col.names=F)
#write.table(spearman_cutoff_indep,file="spearman_cutoff_indep.txt",quote=FALSE,row.names=F,col.names=F)



#saving files for the supplementary file 
supp_data_causal<- causal_data[,c('bcmi_causal','abs_spear_causal','distance')]
supp_data_indep<-indep_data[,c('bcmi_indep','abs_spear_indep','distance')]
write.table(supp_data_causal,"supp_causal_yeast.csv",quote = FALSE,row.names = FALSE)
write.table(supp_data_indep,"supp_indep_yeast.csv",quote = FALSE,row.names = FALSE)
