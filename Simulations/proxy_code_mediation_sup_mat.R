install.packages("gmm")
install.packages("MASS")
install.packages("ivreg")

library(gmm)
library(MASS)
library(ivreg)

rm(list=ls())

source("proxy_fun_red.R")

#Replicate experiments 5-9 in supplementary material.

p.sim <- function(n,nsim,seed,experiment){
  
  results<-matrix(NA,1000,20)
  
  for(i in 1:nsim){
    
    if (experiment==1){
      eps_u<--1;tau_a<--0.3;tau_u<-0;alpha_0u<-0
      l<-dat[,1:2];u<-dat[,3]
      a<-rbinom(n,1,plogis(-l%*%c(0.5,0.5)+alpha_0u*u))
      za_coef=-(1/eps_u)*(tau_a*tau_u+alpha_0u)
      z<-rnorm(n,0.2+za_coef*a+l%*%c(0.2,0.2)+eps_u*u)
      w<-rnorm(n,0.3+l%*%c(0.2,0.2)-0.6*u)
      m<-rnorm(n,tau_a*a-l%*%c(0.5,0.5)+tau_u*u)
      y<-2+2*a+m-l%*%c(1,1)+2*w-u+2*rnorm(n)
    } else if (experiment==2){
      eps_u<--1;tau_a<--0.3;tau_u<-0.4;alpha_0u<--0.4
      l<-dat[,1:2];u<-dat[,3]
      a<-rbinom(n,1,plogis(-l%*%c(0.5,0.5)+alpha_0u*u))
      za_coef=-(1/eps_u)*(tau_a*tau_u+alpha_0u)
      z<-rnorm(n,0.2+za_coef*a+l%*%c(0.2,0.2)+eps_u*u)
      w<-rnorm(n,0.3+l%*%c(0.2,0.2)-0.6*u)
      m<-rnorm(n,tau_a*a-l%*%c(0.5,0.5)+tau_u*u)
      y<-2+2*a+m-l%*%c(1,1)+2*w-u-0.5*z+2*rnorm(n)
    } else if (experiment==3){
      eps_u<--1;tau_a<--0.3;tau_u<-0.4;alpha_0u<--0.4
      l<-dat[,1:2];u<-dat[,3]
      a<-rbinom(n,1,plogis(-l%*%c(0.5,0.5)+alpha_0u*u))
      za_coef=-(1/eps_u)*(tau_a*tau_u+alpha_0u)
      z<-rnorm(n,0.2+za_coef*a+l%*%c(0.2,0.2)+eps_u*u)
      w<-rnorm(n,0.3+l%*%c(0.2,0.2)-0.6*u+0.2*a)
      m<-rnorm(n,tau_a*a-l%*%c(0.5,0.5)+tau_u*u)
      y<-2+2*a+m-l%*%c(1,1)+2*w-u+2*rnorm(n)
    } else if (experiment==4){
      eps_u<--1;tau_a<--0.3;tau_u<-0.4;alpha_0u<--0.4
      l<-dat[,1:2]
      u<-dat[,3]
      a<-rbinom(n,1,plogis(-l%*%c(0.5,0.5)+alpha_0u*u))
      za_coef=-(1/eps_u)*(tau_a*tau_u+alpha_0u)
      z<-rnorm(n,0.2+za_coef*a+l%*%c(0.2,0.2)+eps_u*u)
      w<-rnorm(n,0.3+l%*%c(0.2,0.2)-0.05*u)
      m<-rnorm(n,tau_a*a-l%*%c(0.5,0.5)+tau_u*u)
      y<-2+2*a+m-l%*%c(1,1)+2*w-u+2*rnorm(n)
    } else if (experiment==5){
      eps_u<--0.05;tau_a<--0.1;tau_u<-0.2;alpha_0u<--0.1
      l<-dat[,1:2];u<-dat[,3]
      a<-rbinom(n,1,plogis(-l%*%c(0.5,0.5)+alpha_0u*u))
      za_coef=-(1/eps_u)*(tau_a*tau_u+alpha_0u)
      z<-rnorm(n,0.2+za_coef*a+l%*%c(0.2,0.2)+eps_u*u)
      w<-rnorm(n,0.3+l%*%c(0.2,0.2)-0.6*u)
      m<-rnorm(n,tau_a*a-l%*%c(0.5,0.5)+tau_u*u)
      y<-2+2*a+m-l%*%c(1,1)+2*w-u+2*rnorm(n)
    }
    
    simdata=data.frame(y,a,m,w,z,l)
    

    h1_formula=y~m+l;h0_formula=~l;q1_formula=~m+l;q0_formula=~l
    
    results[i,1]<-coef(lm(y~a+m+l+z+w))[2] 
    results[i,2]<-coef(summary(lm(y~a+m+l+z+w)))[2,2]
    
    p_ipw_fit<-p_ipw_c(simdata,q1_f=q1_formula,q0_f=q0_formula,h0_f=h0_formula,z=z,w=w,out=y,trt=a)
    results[i,3]<-coef(p_ipw_fit)[length(coef(p_ipw_fit))]
    results[i,4]<-sqrt(vcov(p_ipw_fit)[nrow(vcov(p_ipw_fit)),ncol(vcov(p_ipw_fit))])
    
    p_hybrid_fit<-p_hybrid_c(simdata,h1_f=h1_formula,h0_f=h0_formula,q0_f=q0_formula,z=z,w=w,out=y,trt=a)
    results[i,5]<-coef(p_hybrid_fit)[length(coef(p_hybrid_fit))]
    results[i,6]<-sqrt(vcov(p_hybrid_fit)[nrow(vcov(p_hybrid_fit)),ncol(vcov(p_hybrid_fit))])
    
    p_or_fit<-p_or_c(simdata,h1_f=h1_formula,h0_f=h0_formula,q0_f=q0_formula,z=z,w=w,out=y,trt=a)
    results[i,7]<-coef(p_or_fit)[length(coef(p_or_fit))]
    results[i,8]<-sqrt(vcov(p_or_fit)[nrow(vcov(p_or_fit)),ncol(vcov(p_or_fit))])
    
    p_mr_fit<-p_mr_c(simdata,h1_f=h1_formula,h0_f=h0_formula,q1_f=q1_formula,q0_f=q0_formula,z=z,w=w,out=y,trt=a)
    results[i,9]<-coef(p_mr_fit)[length(coef(p_mr_fit))]
    results[i,10]<-sqrt(vcov(p_mr_fit)[nrow(vcov(p_mr_fit)),ncol(vcov(p_mr_fit))])
    
    cat(i,round(c(mean(results[1:i,1]),mean(results[1:i,3]),
                  mean(results[1:i,5]),mean(results[1:i,7]),
                  mean(results[1:i,9])),digits=3),"\n")
  }
  return(results)
}

res5<-p.sim(n=2000,nsim=1000,seed=323,experiment=1)
res6<-p.sim(n=2000,nsim=1000,seed=323,experiment=2) 
res7<-p.sim(n=2000,nsim=1000,seed=323,experiment=3) 
res8<-p.sim(n=2000,nsim=1000,seed=2944,experiment=4)  
res9<-p.sim(n=2000,nsim=1000,seed=2944,experiment=5)  

#Table A1 (supplementary material)

result_mat5<-matrix(NA,5,6)
results_tmp<-res1
seq<-c(1,3,5,7,9)
for(j in seq){
  q<-(j+1)/2
  result_mat5[q,1]<-mean(results_tmp[,j])-2
  result_mat5[q,2]<-median(results_tmp[,j])-2
  result_mat5[q,3]<-mean((results_tmp[,j]-2)^2)
  result_mat5[q,4]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat5[q,5]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat5[q,6]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat6<-matrix(NA,5,6)
results_tmp<-res2
for(j in seq){
  q<-(j+1)/2
  result_mat6[q,1]<-mean(results_tmp[,j])-2
  result_mat6[q,2]<-mean(results_tmp[,j])-2
  result_mat6[q,3]<-mean((results_tmp[,j]-2)^2)
  result_mat6[q,4]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat6[q,5]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat6[q,6]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat7<-matrix(NA,5,6)
results_tmp<-res3
for(j in seq){
  q<-(j+1)/2
  result_mat7[q,1]<-mean(results_tmp[,j])-2
  result_mat7[q,2]<-median(results_tmp[,j])-2
  result_mat7[q,3]<-mean((results_tmp[,j]-2)^2)
  result_mat7[q,4]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat7[q,5]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat7[q,6]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat8<-matrix(NA,5,6)
results_tmp<-res4
for(j in seq){
  q<-(j+1)/2
  result_mat8[q,1]<-mean(results_tmp[,j])-2
  result_mat8[q,2]<-median(results_tmp[,j])-2
  result_mat8[q,3]<-mean((results_tmp[,j]-2)^2)
  result_mat8[q,4]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat8[q,5]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat8[q,6]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat9<-matrix(NA,5,6)
results_tmp<-res5
for(j in seq){
  q<-(j+1)/2
  result_mat9[q,1]<-mean(results_tmp[,j])-2
  result_mat9[q,2]<-median(results_tmp[,j])-2
  result_mat9[q,3]<-mean((results_tmp[,j]-2)^2)
  result_mat9[q,4]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat9[q,5]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat9[q,6]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}

result_mat_full<-round(rbind(result_mat5,result_mat6,result_mat7,result_mat8,result_mat9),digits=2)

estimators<-c("OLS","P-IPW","P-hybrid","P-OR","P-MR","OLS","P-IPW","P-hybrid","P-OR","P-MR","OLS","P-IPW","P-hybrid","P-OR","P-MR","OLS","P-IPW","P-hybrid","P-OR","P-MR","OLS","P-IPW","P-hybrid","P-OR","P-MR")
experiment<-c("5","","","","","6","","","","","7","","","","","8","","","","","9","","","","")
result_mat_fullN<-cbind(experiment,estimators,result_mat_full)
colnames(result_mat_fullN)<-c("Exp","Est","Bias","Med. Bias","MSE","Coverage","Mean Length","Med. Length")
print(xtable(result_mat_fullN),include.rownames=FALSE)



