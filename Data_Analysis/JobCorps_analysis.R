install.packages("mgcv")
install.packages("foreign")
install.packages("MASS")
install.packages("AER")
install.packages("gmm")
install.packages("ivreg")
install.packages("boot")
install.packages("mediation")
install.packages("causalweight")
install.packages("xtable")

library(mgcv)
library(foreign)
library(MASS)
library(AER)
library(gmm)
library(ivreg)
library(boot)
library(mediation)
library(causalweight)
library(xtable)

rm(list=ls())

#Data available from Harvard Dataverse: 
#Huber, Martin, 2020, "JCdata.RData", 
#Replication data for "Direct and indirect effects of continuous treatments based 
#on generalized propensity score weighting", https://doi.org/10.7910/DVN/AJ7Q9B/YZXNTK, 
#Harvard Dataverse, V1

load("JCdata.RData")
source("proxy_fun_complete.R")
attach(JCdata)

#PRE-PROCESSING

x=cbind(female, age, race_white, race_black, race_hispanic,  hgrd, hgrdmissdum, educ_geddiploma, educ_hsdiploma, ntv_engl, marstat_divorced, marstat_separated, marstat_livetogunm, marstat_married, haschldY0, everwkd,  mwearn, hohhd0, peopleathome, peopleathomemissdum, nonres,    g10, g10missdum, g12, g12missdum, hgrd_mum, hgrd_mummissdum, hgrd_dad, hgrd_dadmissdum, work_dad_didnotwork, g2, g5, g7, welfare_child, welfare_childmissdum, h1_fair_poor, h2, h10, h10missdum,  h25, h25missdum, h29, h5, h5missdum, h7, h7missdum, i1, i10, e12, e12missdum, e16, e16missdum, e24usd, e24usdmissdum, e30, e30missdum,  e35, e35missdum,e21,e31, e32, e37, e6_byphone, e8_recruitersoffice, e9ef)
a<-d;a[d>0]<-1

summary(glm(a~x))
summary(lm(y~a+x))

#Check which var is which
summary(e37[d>0]);sd(e37[d>0]);summary(e32[d>0]);sd(e32[d>0])
summary(e21[d>0]);sd(e21[d>0]);summary(e12[d>0&e12missdum==0]);table(e12missdum[d>0])

#Check variables with missing data 
table(hgrdmissdum);table(peopleathomemissdum);table(g10missdum);table(g12missdum)
table(hgrd_mummissdum);table(welfare_childmissdum);table(hgrd_dadmissdum)
table(welfare_childmissdum);table(h10missdum);table(h25missdum);
table(h5missdum);table(h7missdum)

#Bind covariates
l=cbind(female, age, race_white, race_black, race_hispanic,   educ_geddiploma, educ_hsdiploma, ntv_engl, marstat_divorced, marstat_separated, marstat_livetogunm, marstat_married, haschldY0, everwkd,  mwearn, hohhd0, nonres,g10, g10missdum, work_dad_didnotwork, g2, g5, g7, welfare_child, welfare_childmissdum, h1_fair_poor, h2, h29, h5, h5missdum, h7, h7missdum, i1, i10)[e12missdum==0,]

#Create proxies
z<-cbind(e12,e37)[e12missdum==0,]
w<-cbind(e32,e8_recruitersoffice)[e12missdum==0,]

#Create outcome, treatment, mediator
out<-y;y<-out[e12missdum==0]
trt<-a;a<-trt[e12missdum==0]
med<-m;m<-med[e12missdum==0]
df<-data.frame(y,a,m,w,z,l)

#DATA ANALYSIS

#Basic linear model analysis 
out_mod<-lm(y~a+m+l+z+w,data=df)
med_mod<-lm(m~a+l+z+w,data=df)
lin_med<-mediate(med_mod,out_mod,treat = "a", mediator = "m",data=df,boot=FALSE,robustSE=TRUE)
lin_med$z1;lin_med$z0;lin_med$d1;lin_med$d0

#Standard multiply robust approach
proxy_res<-matrix(NA,5,4);se_res<-matrix(NA,5,4)
mr_fit<-mr_c(df,ps1_f=~m+z+w+l,ps0_f=~z+w+l,h1_f=~m+z+w+l,h0_f=~z+w+l,trt=a,out=y)
proxy_res[1,1]<-coef(mr_fit)[length(coef(mr_fit))-3];se_res[1,1]<-sqrt(diag(vcov(mr_fit)))[length(coef(mr_fit))-3]
proxy_res[1,2]<-coef(mr_fit)[length(coef(mr_fit))-2];se_res[1,2]<-sqrt(diag(vcov(mr_fit)))[length(coef(mr_fit))-2]
proxy_res[1,3]<-coef(mr_fit)[length(coef(mr_fit))-1];se_res[1,3]<-sqrt(diag(vcov(mr_fit)))[length(coef(mr_fit))-2]
proxy_res[1,4]<-coef(mr_fit)[length(coef(mr_fit))];se_res[1,4]<-sqrt(diag(vcov(mr_fit)))[length(coef(mr_fit))]

#Proximal approaches
p_ipw_fit<-p_ipw_c(df,q1_f=~m+l,q0_f=~l,h0_f=~l,z=z,w=w,trt=a,out=y)
proxy_res[2,1]<-coef(p_ipw_fit)[length(coef(p_ipw_fit))-3];se_res[2,1]<-sqrt(diag(vcov(p_ipw_fit)))[length(coef(p_ipw_fit))-3]
proxy_res[2,2]<-coef(p_ipw_fit)[length(coef(p_ipw_fit))-2];se_res[2,2]<-sqrt(diag(vcov(p_ipw_fit)))[length(coef(p_ipw_fit))-2]
proxy_res[2,3]<-coef(p_ipw_fit)[length(coef(p_ipw_fit))-1];se_res[2,3]<-sqrt(diag(vcov(p_ipw_fit)))[length(coef(p_ipw_fit))-1]
proxy_res[2,4]<-coef(p_ipw_fit)[length(coef(p_ipw_fit))];se_res[2,4]<-sqrt(diag(vcov(p_ipw_fit)))[length(coef(p_ipw_fit))]

p_hybrid_fit<-p_hybrid_c(df,h1_f=~m+l,q0_f=~l,h0_f=~l,z=z,w=w,trt=a,out=y)
proxy_res[3,1]<-coef(p_hybrid_fit)[length(coef(p_hybrid_fit))-3];se_res[3,1]<-sqrt(diag(vcov(p_hybrid_fit)))[length(coef(p_hybrid_fit))-3]
proxy_res[3,2]<-coef(p_hybrid_fit)[length(coef(p_hybrid_fit))-2];se_res[3,2]<-sqrt(diag(vcov(p_hybrid_fit)))[length(coef(p_hybrid_fit))-2]
proxy_res[3,3]<-coef(p_hybrid_fit)[length(coef(p_hybrid_fit))-1];se_res[3,3]<-sqrt(diag(vcov(p_hybrid_fit)))[length(coef(p_hybrid_fit))-1]
proxy_res[3,4]<-coef(p_hybrid_fit)[length(coef(p_hybrid_fit))];se_res[3,4]<-sqrt(diag(vcov(p_hybrid_fit)))[length(coef(p_hybrid_fit))]

p_or_fit<-p_or_c(df,h1_f=~m,q0_f=~1,h0_f=~1,z=z,w=w,trt=a,out=y)
proxy_res[4,1]<-coef(p_or_fit)[length(coef(p_or_fit))-3];se_res[4,1]<-sqrt(diag(vcov(p_or_fit)))[length(coef(p_or_fit))-3]
proxy_res[4,2]<-coef(p_or_fit)[length(coef(p_or_fit))-2];se_res[4,2]<-sqrt(diag(vcov(p_or_fit)))[length(coef(p_or_fit))-2]
proxy_res[4,3]<-coef(p_or_fit)[length(coef(p_or_fit))-1];se_res[4,3]<-sqrt(diag(vcov(p_or_fit)))[length(coef(p_or_fit))-1]
proxy_res[4,4]<-coef(p_or_fit)[length(coef(p_or_fit))];se_res[4,4]<-sqrt(diag(vcov(p_or_fit)))[length(coef(p_or_fit))]

p_mr_fit<-p_mr_c(df,h1_f=~m+l,q1_f=~m+l,q0_f=~l,h0_f=~l,z=z,w=w,trt=a,out=y)
proxy_res[5,1]<-coef(p_mr_fit)[length(coef(p_mr_fit))-3];se_res[5,1]<-sqrt(diag(vcov(p_mr_fit)))[length(coef(p_mr_fit))-3]
proxy_res[5,2]<-coef(p_mr_fit)[length(coef(p_mr_fit))-2];se_res[5,2]<-sqrt(diag(vcov(p_mr_fit)))[length(coef(p_mr_fit))-2]
proxy_res[5,3]<-coef(p_mr_fit)[length(coef(p_mr_fit))-1];se_res[5,3]<-sqrt(diag(vcov(p_mr_fit)))[length(coef(p_mr_fit))-1]
proxy_res[5,4]<-coef(p_mr_fit)[length(coef(p_mr_fit))];se_res[5,4]<-sqrt(diag(vcov(p_mr_fit)))[length(coef(p_mr_fit))]

#Table 2 (main manuscript)

lat_tab<-matrix(NA,4,15)
for(j in 1:4){
lat_tab[j,1]<-proxy_res[1,j]
lat_tab[j,2]<-proxy_res[1,j]-1.96*se_res[1,j]
lat_tab[j,3]<-proxy_res[1,j]+1.96*se_res[1,j]
lat_tab[j,4]<-proxy_res[2,j]
lat_tab[j,5]<-proxy_res[2,j]-1.96*se_res[2,j]
lat_tab[j,6]<-proxy_res[2,j]+1.96*se_res[2,j]
lat_tab[j,7]<-proxy_res[3,j]
lat_tab[j,8]<-proxy_res[3,j]-1.96*se_res[3,j]
lat_tab[j,9]<-proxy_res[3,j]+1.96*se_res[3,j]
lat_tab[j,10]<-proxy_res[4,j]
lat_tab[j,11]<-proxy_res[4,j]-1.96*se_res[4,j]
lat_tab[j,12]<-proxy_res[4,j]+1.96*se_res[4,j]
lat_tab[j,13]<-proxy_res[5,j]
lat_tab[j,14]<-proxy_res[5,j]-1.96*se_res[5,j]
lat_tab[j,15]<-proxy_res[5,j]+1.96*se_res[5,j]
}
rownames(lat_tab)=c("Direct Effect 0","Indirect Effect 1","Direct Effect 1","Indirect Effect 0")
colnames(lat_tab)=c("S-MR","CI LL","CI UL","P-IPW","CI LL","CI UL","P-Hybrid","CI LL","CI UL","P-OR","CI LL","CI UL","P-MR","CI LL","CI UL")
print(xtable(lat_tab,digits=4), sanitize.colnames.function = function(x) {x})

lat_tab_red<-lat_tab[,c(1:3,13:15)]
print(xtable(lat_tab_red,digits=4), sanitize.colnames.function = function(x) {x})

#Total effects
#Standard approach
coef(mr_fit)[length(coef(mr_fit))-4]
coef(mr_fit)[length(coef(mr_fit))-4]-sqrt(diag(vcov(mr_fit)))[length(coef(mr_fit))-4]
coef(mr_fit)[length(coef(mr_fit))-4]+sqrt(diag(vcov(mr_fit)))[length(coef(mr_fit))-4]

#Proximal approach
coef(p_mr_fit)[length(coef(p_mr_fit))-4]
coef(p_mr_fit)[length(coef(p_mr_fit))-4]-sqrt(diag(vcov(p_mr_fit)))[length(coef(p_mr_fit))-4]
coef(p_mr_fit)[length(coef(p_mr_fit))-4]+sqrt(diag(vcov(p_mr_fit)))[length(coef(p_mr_fit))-4]
