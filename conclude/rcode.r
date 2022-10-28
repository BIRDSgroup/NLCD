library(cit)
L=rbinom(1000,1,0.5)
A=6*L + rnorm(1000,0,1)
B[L==1]=-3*(A[L==1]-6)*(A[L==1]-6)+rnorm(length(A[L==1]),0,1)+10
B[L==0]=3*(A[L==0])*(A[L==0])+rnorm(length(A[L==0]),0,1)-24 #this is wrong
lm(B~A+L)
points(A,Bstar1,color='red')
fit=lm(Bstar~A+L)
Bstar[L==1]=sample(B[L==1])
ggplot(df,aes(A,B,color=L))+geom_point()
B0=predict(fit,input0_actual)
input1_actual[L==0,2]=1
input=data.frame(A,L)