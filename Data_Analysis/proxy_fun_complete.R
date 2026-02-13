q0_fun<-function(gamma,int.cov_z,a){
  return(1+exp(int.cov_z%*%gamma))
}
moments_q0t<-function(gamma,data){
  int.cov_z<-data.matrix(cbind(data$cov_q0,data$z));int.cov_w<-data.matrix(cbind(data$cov_q0,data$w))
  a<-cbind(data$a)
  ee<-int.cov_w*as.vector(a*q0_fun(gamma,int.cov_z=int.cov_z,a=a)-1)
  return(cbind(ee))
}
moments_q0c<-function(gamma,data){
  int.cov_z<-data.matrix(cbind(data$cov_q0,data$z));int.cov_w<-data.matrix(cbind(data$cov_q0,data$w))
  a<-cbind(data$a)
  ee<-int.cov_w*as.vector((1-a)*q0_fun(gamma,int.cov_z=int.cov_z,a=a)-1)
  return(cbind(ee))
}

moments_q1_t<-function(gamma,data){
  int.cov_z<-data.matrix(cbind(data$cov_q1,data$z))
  int.cov_w<-data.matrix(cbind(data$cov_q1,data$w))
  a<-cbind(data$a);q0c<-as.vector(cbind(data$q0c))
  ee<-int.cov_w*as.vector(a*(exp(int.cov_z%*%gamma))-(1-a)*q0c)
  return(cbind(ee))
}

moments_q1_c<-function(gamma,data){
  int.cov_z<-data.matrix(cbind(data$cov_q1,data$z))
  int.cov_w<-data.matrix(cbind(data$cov_q1,data$w))
  a<-cbind(data$a);q0t<-as.vector(cbind(data$q0t))
  ee<-int.cov_w*as.vector((1-a)*q0t*(exp(int.cov_z%*%gamma))-a*q0t)
  return(cbind(ee))
}

p_or_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_h1<-data.matrix(cbind(data$cov_h1));cov_h0<-data.matrix(cbind(data$cov_h0))
  cov_q0<-data.matrix(cbind(data$cov_q0))

  h1t_pos <- 1:ncol(cbind(cov_h1,w))
  h1c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h1,w))))
  h0t_pos <- (max(h1c_pos) + 1):(max(h1c_pos) + (ncol(cbind(cov_h0,w))))
  h0c_pos <- (max(h0t_pos) + 1):(max(h0t_pos) + (ncol(cbind(cov_h0,w))))
  h2t_pos <- (max(h0c_pos) + 1):(max(h0c_pos) + (ncol(cbind(cov_h0,w))))
  h2c_pos <- (max(h2t_pos) + 1):(max(h2t_pos) + (ncol(cbind(cov_h0,w))))
  q0t_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_q0,z))))
  q0c_pos <- (max(q0t_pos) + 1):(max(q0t_pos) + (ncol(cbind(cov_q0,z))))
  te_pos <- (max(q0c_pos) + 1)

  h1t<-cbind(cov_h1,w)%*%psi[h1t_pos];h1c<-cbind(cov_h1,w)%*%psi[h1c_pos]
  h0t<-cbind(cov_h0,w)%*%psi[h0t_pos];h0c<-cbind(cov_h0,w)%*%psi[h0c_pos]
  h2t<-cbind(cov_h0,w)%*%psi[h2t_pos];h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  q0t <- as.vector(q0_fun(gamma=psi[q0t_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))

  m1<-cbind(cov_h1,z)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h1,z)*as.vector((1-a)*(y-h1c))
  m3<-cbind(cov_h0,z)*as.vector(a*(h1c-h0t))
  m4<-cbind(cov_h0,z)*as.vector((1-a)*(h1t-h0c))
  m5<-cbind(cov_h0,z)*as.vector(a*(y-h2t))
  m6<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m7<-cbind(cov_q0,w)*as.vector(a*q0t-1)
  m8<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m9<-h0c-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] #E{Y(1,M(0))-Y(0)}
  m10<-(a*q0t*(y-h2t)+h2t)-h0t-psi[te_pos+1]  #E{Y(1)-Y(0,M(1))}
  m11<-(a*q0t*(y-h2t)+h2t)-h0c-psi[te_pos+2]  #E{Y(1)-Y(1,M(0))}
  m12<-h0t-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos+3]  #E{Y(0,M(1))-Y(0)}
  return(cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12))
}

p_or_c<-function(data,h1_f,h0_f,q0_f,w,z,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_q0<-model.matrix(q0_f,data=data);y=out;a=trt

  h1t_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==1)
  h1t <- data.matrix(cbind(cov_h1,w))%*%h1t_fit$coefficients
  h1c_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==0)
  h1c <- data.matrix(cbind(cov_h1,w))%*%h1c_fit$coefficients
  
  h0c_fit<-ivreg(h1t~-1+cov_h0+w|cov_h0+z,subset=a==0)
  h0t_fit<-ivreg(h1c~-1+cov_h0+w|cov_h0+z,subset=a==1)

  h2t_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==1)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)

  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0t_fit<-gmm(moments_q0t, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  
  init <- c(coef(h1t_fit),coef(h1c_fit),coef(h0t_fit),coef(h0c_fit),
            coef(h2t_fit),coef(h2c_fit),coef(q0t_fit),coef(q0c_fit),
            rep(0,4));nparms=length(init)
            
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_h1=cov_h1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_or_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-4)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

#p_or_c(simdata,h1_f=~m+l,h0_f=~l,q0_f=~l,z=z,w=w,trt=a,out=y)


p_hybrid_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_h1<-data.matrix(cbind(data$cov_h1));cov_h0<-data.matrix(cbind(data$cov_h0))
  cov_q0<-data.matrix(cbind(data$cov_q0))
  
  h1t_pos <- 1:ncol(cbind(cov_h1,w))
  h1c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h1,w))))
  h2t_pos <- (max(h1c_pos) + 1):(max(h1c_pos) + (ncol(cbind(cov_h0,w))))
  h2c_pos <- (max(h2t_pos) + 1):(max(h2t_pos) + (ncol(cbind(cov_h0,w))))
  q0t_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_q0,z))))
  q0c_pos <- (max(q0t_pos) + 1):(max(q0t_pos) + (ncol(cbind(cov_q0,z))))
  te_pos <- (max(q0c_pos) + 1)
  
  h1t<-cbind(cov_h1,w)%*%psi[h1t_pos];h1c<-cbind(cov_h1,w)%*%psi[h1c_pos]
  h2t<-cbind(cov_h0,w)%*%psi[h2t_pos];h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  q0t <- as.vector(q0_fun(gamma=psi[q0t_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  
  m1<-cbind(cov_h1,z)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h1,z)*as.vector((1-a)*(y-h1c))
  m3<-cbind(cov_h0,z)*as.vector(a*(y-h2t))
  m4<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m5<-cbind(cov_q0,w)*as.vector(a*q0t-1)
  m6<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m7<-(1-a)*q0c*h1t-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] #E{Y(1,M(0))-Y(0)}
  m8<-(a*q0t*(y-h2t)+h2t)-a*q0t*h1c-psi[te_pos+1]  #E{Y(1)-Y(0,M(1))}
  m9<-(a*q0t*(y-h2t)+h2t)-(1-a)*q0c*h1t-psi[te_pos+2]  #E{Y(1)-Y(1,M(0))}
  m10<-a*q0t*h1c-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos+3]  #E{Y(0,M(1))-Y(0)}
  return(cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10))
}

p_hybrid_c<-function(data,h1_f,h0_f,q0_f,w,z,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_q0<-model.matrix(q0_f,data=data);y=out;a=trt
  
  h1t_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==1)
  h1c_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==0)
  h2t_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==1)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0t_fit<-gmm(moments_q0t, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  
  init <- c(coef(h1t_fit),coef(h1c_fit),coef(h2t_fit),coef(h2c_fit),
            coef(q0t_fit),coef(q0c_fit),rep(0,4));nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_h1=cov_h1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_hybrid_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-4)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

#p_hybrid_c(simdata,h1_f=~m+l,h0_f=~l,q0_f=~l,z=z,w=w,out=y,trt=a)

p_ipw_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_q1<-data.matrix(cbind(data$cov_q1));cov_h0<-data.matrix(cbind(data$cov_h0))
  cov_q0<-data.matrix(cbind(data$cov_q0))
  
  q0t_pos <- 1:ncol(cbind(cov_q0,z))
  q0c_pos <- (max(q0t_pos) + 1):(max(q0t_pos) + (ncol(cbind(cov_q0,z))))
  q1t_pos <- (max(q0c_pos) + 1):(max(q0c_pos) + (ncol(cbind(cov_q1,z))))
  q1c_pos <- (max(q1t_pos) + 1):(max(q1t_pos) + (ncol(cbind(cov_q1,z))))
  h2t_pos <- (max(q1c_pos) + 1):(max(q1c_pos) + (ncol(cbind(cov_h0,w))))
  h2c_pos <- (max(h2t_pos) + 1):(max(h2t_pos) + (ncol(cbind(cov_h0,w))))
  te_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + 1)
  
  q0t <- as.vector(q0_fun(gamma=psi[q0t_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q1t <- as.vector(exp(cbind(cov_q1,z)%*%psi[q1t_pos]));q1c<-as.vector(exp(cbind(cov_q1,z)%*%psi[q1c_pos]))
  h2t<-cbind(cov_h0,w)%*%psi[h2t_pos];h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]

  m1<-cbind(cov_q0,w)*as.vector(a*q0t-1)
  m2<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m3<-cbind(cov_q1,w)*(a*q0c*q1t-(1-a)*q0c)
  m4<-cbind(cov_q1,w)*((1-a)*q0t*q1c-a*q0t)
  m5<-cbind(cov_h0,z)*as.vector(a*(y-h2t))
  m6<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m7<-a*q0c*q1t*y-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] #E{Y(1,M(0))-Y(0)}
  m8<-(a*q0t*(y-h2t)+h2t)-(1-a)*q0t*q1c*y-psi[te_pos+1]  #E{Y(1)-Y(0,M(1))}
  m9<-(a*q0t*(y-h2t)+h2t)-a*q0c*q1t*y-psi[te_pos+2]  #E{Y(1)-Y(1,M(0))}
  m10<-(1-a)*q0t*q1c*y-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos+3]  #E{Y(0,M(1))-Y(0)}
  return(cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10))
}

p_ipw_c<-function(data,q1_f,q0_f,h0_f,w,z,trt,out){
  cov_q1<-model.matrix(q1_f,data=data);cov_q0<-model.matrix(q0_f,data=data)
  cov_h0<-model.matrix(h0_f,data=data);y=out;a=trt
  
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0t_fit<-gmm(moments_q0t, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0t<-as.vector(q0_fun(gamma=coef(q0t_fit),int.cov_z=cbind(cov_q0,z),a=a))
  q0c<-as.vector(q0_fun(gamma=coef(q0c_fit),int.cov_z=cbind(cov_q0,z),a=a))
  q1_data<-list(cov_q1=cov_q1,a=a,z=z,w=w,q0t=q0t,q0c=q0c);init<-rep(0,dim(cov_q1)[2]+dim(as.matrix(w))[2])
  q1t_fit<-gmm(moments_q1_t, x = q1_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q1c_fit<-gmm(moments_q1_c, x = q1_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  h2t_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==1)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  
  init <- c(coef(q0t_fit),coef(q0c_fit),coef(q1t_fit),coef(q1c_fit),coef(h2t_fit),coef(h2c_fit),rep(0,4));nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_q1=cov_q1,cov_h0=cov_h0,cov_q0=cov_q0)

  gmm(p_ipw_fun_c, x = gmm_data, t0 = init,
       eqConst=c(1:(nparms-4)),eqConstFullVcov=TRUE,
       crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

#p_ipw_c(simdata,h0_f=~l,q0_f=~l,q1_f=~m+l,z=z,w=w,out=y,trt=a)


p_mr_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_h1<-data.matrix(cbind(data$cov_h1));cov_q1<-data.matrix(cbind(data$cov_q1))
  cov_h0<-data.matrix(cbind(data$cov_h0));cov_q0<-data.matrix(cbind(data$cov_q0))
  
  h1t_pos <- 1:ncol(cbind(cov_h1,w))
  h1c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h1,w))))
  h0t_pos <- (max(h1c_pos) + 1):(max(h1c_pos) + (ncol(cbind(cov_h0,w))))
  h0c_pos <- (max(h0t_pos) + 1):(max(h0t_pos) + (ncol(cbind(cov_h0,w))))
  h2t_pos <- (max(h0c_pos) + 1):(max(h0c_pos) + (ncol(cbind(cov_h0,w))))
  h2c_pos <- (max(h2t_pos) + 1):(max(h2t_pos) + (ncol(cbind(cov_h0,w))))
  q0t_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_q0,z))))
  q0c_pos <- (max(q0t_pos) + 1):(max(q0t_pos) + (ncol(cbind(cov_q0,z))))
  q1t_pos <- (max(q0c_pos) + 1):(max(q0c_pos) + (ncol(cbind(cov_q1,z))))
  q1c_pos <- (max(q1t_pos) + 1):(max(q1t_pos) + (ncol(cbind(cov_q1,z))))
  te_pos <- (max(q1c_pos) + 1)
  
  h1t<-cbind(cov_h1,w)%*%psi[h1t_pos];h1c<-cbind(cov_h1,w)%*%psi[h1c_pos]
  h0t<-cbind(cov_h0,w)%*%psi[h0t_pos];h0c<-cbind(cov_h0,w)%*%psi[h0c_pos]
  h2t<-cbind(cov_h0,w)%*%psi[h2t_pos];h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  q0t <- as.vector(q0_fun(gamma=psi[q0t_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q1t <- as.vector(exp(cbind(cov_q1,z)%*%psi[q1t_pos]));q1c<-as.vector(exp(cbind(cov_q1,z)%*%psi[q1c_pos]))
  
  m1<-cbind(cov_h1,z)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h1,z)*as.vector((1-a)*(y-h1c))
  m3<-cbind(cov_h0,z)*as.vector(a*(h1c-h0t))
  m4<-cbind(cov_h0,z)*as.vector((1-a)*(h1t-h0c))
  m5<-cbind(cov_h0,z)*as.vector(a*(y-h2t))
  m6<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m7<-cbind(cov_q0,w)*as.vector(a*q0t-1)
  m8<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m9<-cbind(cov_q1,w)*(a*q0c*q1t-(1-a)*q0c)
  m10<-cbind(cov_q1,w)*((1-a)*q0t*q1c-a*q0t)
  m11<-(a*q0t*(y-h2t)+h2t)-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] 
  m12<-(a*q0c*q1t*(y-h1t)+(1-a)*q0c*(h1t-h0c)+h0c)-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos+1] 
  m13<-(a*q0t*(y-h2t)+h2t)-((1-a)*q0t*q1c*(y-h1c)+a*q0t*(h1c-h0t)+h0t)-psi[te_pos+2]
  m14<-(a*q0t*(y-h2t)+h2t)-(a*q0c*q1t*(y-h1t)+(1-a)*q0c*(h1t-h0c)+h0c)-psi[te_pos+3]
  m15<-((1-a)*q0t*q1c*(y-h1c)+a*q0t*(h1c-h0t)+h0t)-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos+4]
  return(cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15))
}

p_mr_c<-function(data,h1_f,h0_f,q1_f,q0_f,w,z,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_q1<-model.matrix(q1_f,data=data);cov_q0<-model.matrix(q0_f,data=data)
  y=out;a=trt
  
  h1t_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==1)
  h1t <- data.matrix(cbind(cov_h1,w))%*%h1t_fit$coefficients
  h1c_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==0)
  h1c <- data.matrix(cbind(cov_h1,w))%*%h1c_fit$coefficients
  h0c_fit<-ivreg(h1t~-1+cov_h0+w|cov_h0+z,subset=a==0)
  h0t_fit<-ivreg(h1c~-1+cov_h0+w|cov_h0+z,subset=a==1)
  h2t_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==1)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0t_fit<-gmm(moments_q0t, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0t<-as.vector(q0_fun(gamma=coef(q0t_fit),int.cov_z=cbind(cov_q0,z),a=a))
  q0c<-as.vector(q0_fun(gamma=coef(q0c_fit),int.cov_z=cbind(cov_q0,z),a=a))
  q1_data<-list(cov_q1=cov_q1,a=a,z=z,w=w,q0t=q0t,q0c=q0c);init<-rep(0,dim(cov_q1)[2]+dim(as.matrix(w))[2])
  q1t_fit<-gmm(moments_q1_t, x = q1_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q1c_fit<-gmm(moments_q1_c, x = q1_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  
  init <- c(coef(h1t_fit),coef(h1c_fit),coef(h0t_fit),coef(h0c_fit),
            coef(h2t_fit),coef(h2c_fit),coef(q0t_fit),coef(q0c_fit),
            coef(q1t_fit),coef(q1c_fit),rep(0,5));nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_h1=cov_h1,cov_q1=cov_q1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_mr_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-5)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

mr_fun_c<-function(psi,data){
  cov_h1<-data.matrix(cbind(data$cov_h1));cov_ps1<-data.matrix(cbind(data$cov_ps1))
  cov_h0<-data.matrix(cbind(data$cov_h0));cov_ps0<-data.matrix(cbind(data$cov_ps0))
  
  h1t_pos <- 1:ncol(cbind(cov_h1))
  h1c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h1))))
  h0t_pos <- (max(h1c_pos) + 1):(max(h1c_pos) + (ncol(cbind(cov_h0))))
  h0c_pos <- (max(h0t_pos) + 1):(max(h0t_pos) + (ncol(cbind(cov_h0))))
  h2t_pos <- (max(h0c_pos) + 1):(max(h0c_pos) + (ncol(cbind(cov_h0))))
  h2c_pos <- (max(h2t_pos) + 1):(max(h2t_pos) + (ncol(cbind(cov_h0))))
  ps0_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_ps0))))
  ps1_pos <- (max(ps0_pos) + 1):(max(ps0_pos) + (ncol(cbind(cov_ps1))))
  te_pos <- (max(ps1_pos) + 1)
  
  h1t<-cbind(cov_h1)%*%psi[h1t_pos];h1c<-cbind(cov_h1)%*%psi[h1c_pos]
  h0t<-cbind(cov_h0)%*%psi[h0t_pos];h0c<-cbind(cov_h0)%*%psi[h0c_pos]
  h2t<-cbind(cov_h0)%*%psi[h2t_pos];h2c<-cbind(cov_h0)%*%psi[h2c_pos]
  ps0<-plogis(cov_ps0%*%psi[ps0_pos]);ps1<-plogis(cov_ps1%*%psi[ps1_pos])
  w_t<-((1-ps1)/((1-ps0)*ps1));w_c<-(ps1/(ps0*(1-ps1)))
  dr1<-(a/ps0)*(y-h2t)+h2t;dr0<-((1-a)/(1-ps0))*(y-h2c)+h2c
  
  m1<-cbind(cov_h1)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h1)*as.vector((1-a)*(y-h1c))
  m3<-cbind(cov_h0)*as.vector(a*(h1c-h0t))
  m4<-cbind(cov_h0)*as.vector((1-a)*(h1t-h0c))
  m5<-cbind(cov_h0)*as.vector(a*(y-h2t))
  m6<-cbind(cov_h0)*as.vector((1-a)*(y-h2c))
  m7<-cbind(cov_ps0)*as.vector(a-ps0)
  m8<-cbind(cov_ps1)*as.vector(a-ps1)
  m9<-dr1-dr0-psi[te_pos] 
  m10<-(a*w_t*(y-h1t)+((1-a)/(1-ps0))*(h1t-h0c)+h0c)-dr0-psi[te_pos+1] 
  m11<-dr1-((1-a)*w_c*(y-h1c)+(a/ps0)*(h1c-h0t)+h0t)-psi[te_pos+2] 
  m12<-dr1-(a*w_t*(y-h1t)+((1-a)/(1-ps0))*(h1t-h0c)+h0c)-psi[te_pos+3] 
  m13<-((1-a)*w_c*(y-h1c)+(a/ps0)*(h1c-h0t)+h0t)-dr0-psi[te_pos+4] 
  return(cbind(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13))
}

mr_c<-function(data,h1_f,h0_f,ps1_f,ps0_f,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_ps1<-model.matrix(ps1_f,data=data);cov_ps0<-model.matrix(ps0_f,data=data)
  y=out;a=trt

  h1t_fit<-lm(y~-1+cov_h1,subset=a==1)
  h1t <- data.matrix(cbind(cov_h1))%*%h1t_fit$coefficients
  h1c_fit<-lm(y~-1+cov_h1,subset=a==0)
  h1c <- data.matrix(cbind(cov_h1))%*%h1c_fit$coefficients
  h0c_fit<-lm(h1t~-1+cov_h0,subset=a==0);h0t_fit<-lm(h1c~-1+cov_h0,subset=a==1)
  h2t_fit<-lm(y~-1+cov_h0,subset=a==1);h2c_fit<-lm(y~-1+cov_h0,subset=a==0)
  ps0_fit<-glm(a~-1+cov_ps0,family="binomial")
  ps1_fit<-glm(a~-1+cov_ps1,family="binomial")
  
  init <- c(coef(h1t_fit),coef(h1c_fit),coef(h0t_fit),coef(h0c_fit),
            coef(h2t_fit),coef(h2c_fit),coef(ps0_fit),coef(ps1_fit),
            rep(0,5));nparms=length(init)
  
  gmm_data<-list(y=y,a=a,cov_h1=cov_h1,cov_ps1=cov_ps1,cov_h0=cov_h0,cov_ps0=cov_ps0)
  
  gmm(mr_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-5)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}