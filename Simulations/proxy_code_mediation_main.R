install.packages("gmm")
install.packages("MASS")
install.packages("ivreg")

library(gmm)
library(MASS)
library(ivreg)

rm(list=ls())

source("proxy_fun_red.R")

#Replicate experiments 1-4 in main manuscript.

p.sim <- function(n,nsim,seed,experiment){
  
  results<-matrix(NA,1000,10)
  
  for(i in 1:nsim){
    
    set.seed(i*seed)
    sigma<-matrix(c(0.5^2,0,0.05,0,0.5^2,0.05,0.05,0.05,1),nrow=3,ncol=3)
    dat<-mvrnorm(n,c(0.25,0.25,0),sigma)
    l<-dat[,1:2]
    u<-dat[,3]
    a<-rbinom(n,1,plogis(-l%*%c(0.5,0.5)-0.4*u))
    za_coef=-0.3*0.4-0.4
    z<-rnorm(n,0.2+za_coef*a+l%*%c(0.2,0.2)-u) #Make z dependent on A.
    w<-rnorm(n,0.3+l%*%c(0.2,0.2)-0.6*u)
    m<-rnorm(n,-0.3*a-l%*%c(0.5,0.5)+0.4*u)
    y<-2+2*a+m-l%*%c(1,1)+2*w-u+2*rnorm(n)
    cor(z,w)
    
    simdata=data.frame(y,a,m,w,z,l)
    
    l_mis<-cbind(abs(l[,1])^(0.5),abs(l[,2])^(0.5))
    
    if (experiment==1){ #All models correct: M_int
      h1_formula=y~m+l;h0_formula=~l;q1_formula=~m+l;q0_formula=~l
    } else if (experiment==2) { #Both Q models wrong
      h1_formula=y~m+l;h0_formula=~l;q1_formula=~m+l_mis;q0_formula=~l_mis
    } else if (experiment==3) { #H0 and Q1 wrong
      h1_formula=y~m+l;h0_formula=~l_mis;q1_formula=~m+l_mis;q0_formula=~l
    } else if (experiment==4) { #Both H models wrong
      h1_formula=y~m+l_mis;h0_formula=~l_mis;q1_formula=~m+l;q0_formula=~l
    } 
    
    results[i,1]<-coef(lm(y~a+m+l+z+w))[2] #Naive direct effect estimate
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

res1<-p.sim(n=2000,nsim=1000,seed=323,experiment=1)
res2<-p.sim(n=2000,nsim=1000,seed=254,experiment=2)
res3<-p.sim(n=2000,nsim=1000,seed=499,experiment=3)
res4<-p.sim(n=2000,nsim=1000,seed=7,experiment=4)

#Table 1 (main manuscript)

result_mat1<-matrix(NA,4,5)
results_tmp<-res1
seq<-c(3,5,7,9)
for(j in seq){
  q<-(j+1)/2-1
  result_mat1[q,1]<-mean(results_tmp[,j])-2
  result_mat1[q,2]<-mean((results_tmp[,j]-2)^2)
  result_mat1[q,3]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat1[q,4]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat1[q,5]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat2<-matrix(NA,4,5)
results_tmp<-res2
for(j in seq){
  q<-(j+1)/2-1
  result_mat2[q,1]<-mean(results_tmp[,j])-2
  result_mat2[q,2]<-mean((results_tmp[,j]-2)^2)
  result_mat2[q,3]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat2[q,4]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat2[q,5]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat3<-matrix(NA,4,5)
results_tmp<-res3
for(j in seq){
  q<-(j+1)/2-1
  result_mat3[q,1]<-mean(results_tmp[,j])-2
  result_mat3[q,2]<-mean((results_tmp[,j]-2)^2)
  result_mat3[q,3]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat3[q,4]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat3[q,5]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat4<-matrix(NA,4,5)
results_tmp<-res4
for(j in seq){
  q<-(j+1)/2-1
  result_mat4[q,1]<-mean(results_tmp[,j])-2
  result_mat4[q,2]<-mean((results_tmp[,j]-2)^2)
  result_mat4[q,3]<-mean(((results_tmp[,j]-1.96*results_tmp[,j+1])<2)&((results_tmp[,j]+1.96*results_tmp[,j+1])>2))
  result_mat4[q,4]<-mean((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
  result_mat4[q,5]<-median((results_tmp[,j]+1.96*results_tmp[,j+1])-(results_tmp[,j]-1.96*results_tmp[,j+1]))
}
result_mat_full<-round(rbind(result_mat1,result_mat2,result_mat3,result_mat4),digits=2)

estimators<-c("P-IPW","P-hybrid","P-OR","P-MR","P-IPW","P-hybrid","P-OR","P-MR","P-IPW","P-hybrid","P-OR","P-MR","P-IPW","P-hybrid","P-OR","P-MR")
experiment<-c("1","","","","2","","","","3","","","","4","","","")
result_mat_fullN<-cbind(experiment,estimators,result_mat_full)
colnames(result_mat_fullN)<-c("Exp","Est","Bias","MSE","Coverage","Mean Length","Med. Length")
print(xtable(result_mat_fullN),include.rownames=FALSE)


