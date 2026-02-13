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
  ee<-int.cov_w*as.vector(a*q0c*(exp(int.cov_z%*%gamma))-(1-a)*q0c)
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
  h0c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h0,w))))
  h2c_pos <- (max(h0c_pos) + 1):(max(h0c_pos) + (ncol(cbind(cov_h0,w))))
  q0c_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_q0,z))))
  te_pos <- (max(q0c_pos) + 1)
  
  h1t<-cbind(cov_h1,w)%*%psi[h1t_pos]
  h0c<-cbind(cov_h0,w)%*%psi[h0c_pos]
  h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  
  m1<-cbind(cov_h1,z)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h0,z)*as.vector((1-a)*(h1t-h0c))
  m3<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m4<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m5<-h0c-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] #E{Y(1,M(0))-Y(0)}
  return(cbind(m1,m2,m3,m4,m5))
}

p_or_c<-function(data,h1_f,h0_f,q0_f,w,z,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_q0<-model.matrix(q0_f,data=data);y=out;a=trt
  
  h1t_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==1)
  h1t <- data.matrix(cbind(cov_h1,w))%*%h1t_fit$coefficients
  h0c_fit<-ivreg(h1t~-1+cov_h0+w|cov_h0+z,subset=a==0)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  
  init <- c(coef(h1t_fit),coef(h0c_fit),coef(h2c_fit),coef(q0c_fit),0);nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_h1=cov_h1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_or_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-1)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

p_hybrid_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_h1<-data.matrix(cbind(data$cov_h1));cov_h0<-data.matrix(cbind(data$cov_h0))
  cov_q0<-data.matrix(cbind(data$cov_q0))
  
  h1t_pos <- 1:ncol(cbind(cov_h1,w))
  h2c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h0,w))))
  q0c_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_q0,z))))
  te_pos <- (max(q0c_pos) + 1)
  
  h1t<-cbind(cov_h1,w)%*%psi[h1t_pos];h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  
  m1<-cbind(cov_h1,z)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m3<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m4<-(1-a)*q0c*h1t-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] #E{Y(1,M(0))-Y(0)}
  return(cbind(m1,m2,m3,m4))
}

p_hybrid_c<-function(data,h1_f,h0_f,q0_f,w,z,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_q0<-model.matrix(q0_f,data=data);y=out;a=trt
  
  h1t_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==1)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  
  init <- c(coef(h1t_fit),coef(h2c_fit),coef(q0c_fit),0);nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_h1=cov_h1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_hybrid_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-1)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

p_ipw_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_q1<-data.matrix(cbind(data$cov_q1));cov_h0<-data.matrix(cbind(data$cov_h0))
  cov_q0<-data.matrix(cbind(data$cov_q0))
  
  q0c_pos <- 1:ncol(cbind(cov_q0,z))
  q1t_pos <- (max(q0c_pos) + 1):(max(q0c_pos) + (ncol(cbind(cov_q1,z))))
  h2c_pos <- (max(q1t_pos) + 1):(max(q1t_pos) + (ncol(cbind(cov_h0,w))))
  te_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + 1)
  
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q1t <- as.vector(exp(cbind(cov_q1,z)%*%psi[q1t_pos]))
  h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  
  m1<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m2<-cbind(cov_q1,w)*(a*q0c*q1t-(1-a)*q0c)
  m3<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m4<-a*q0c*q1t*y-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] #E{Y(1,M(0))-Y(0)}
  return(cbind(m1,m2,m3,m4))
}

p_ipw_c<-function(data,q1_f,q0_f,h0_f,w,z,trt,out){
  cov_q1<-model.matrix(q1_f,data=data);cov_q0<-model.matrix(q0_f,data=data)
  cov_h0<-model.matrix(h0_f,data=data);y=out;a=trt
  
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0c<-as.vector(q0_fun(gamma=coef(q0c_fit),int.cov_z=cbind(cov_q0,z),a=a))
  q1_data<-list(cov_q1=cov_q1,a=a,z=z,w=w,q0c=q0c);init<-rep(0,dim(cov_q1)[2]+dim(as.matrix(w))[2])
  q1t_fit<-gmm(moments_q1_t, x = q1_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  
  init <- c(coef(q0c_fit),coef(q1t_fit),coef(h2c_fit),0);nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_q1=cov_q1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_ipw_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-1)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

p_mr_fun_c<-function(psi,data){
  y=as.numeric(data$y);a=as.numeric(data$a)
  z=data.matrix(cbind(data$z));w=data.matrix(cbind(data$w))
  cov_h1<-data.matrix(cbind(data$cov_h1));cov_q1<-data.matrix(cbind(data$cov_q1))
  cov_h0<-data.matrix(cbind(data$cov_h0));cov_q0<-data.matrix(cbind(data$cov_q0))
  
  h1t_pos <- 1:ncol(cbind(cov_h1,w))
  h0c_pos <- (max(h1t_pos) + 1):(max(h1t_pos) + (ncol(cbind(cov_h0,w))))
  h2c_pos <- (max(h0c_pos) + 1):(max(h0c_pos) + (ncol(cbind(cov_h0,w))))
  q0c_pos <- (max(h2c_pos) + 1):(max(h2c_pos) + (ncol(cbind(cov_q0,z))))
  q1t_pos <- (max(q0c_pos) + 1):(max(q0c_pos) + (ncol(cbind(cov_q1,z))))
  te_pos <- (max(q1t_pos) + 1)
  
  h1t<-cbind(cov_h1,w)%*%psi[h1t_pos];h0c<-cbind(cov_h0,w)%*%psi[h0c_pos]
  h2c<-cbind(cov_h0,w)%*%psi[h2c_pos]
  q0c <- as.vector(q0_fun(gamma=psi[q0c_pos],int.cov_z=cbind(cov_q0,z),a=a))
  q1t <- as.vector(exp(cbind(cov_q1,z)%*%psi[q1t_pos]))
  
  m1<-cbind(cov_h1,z)*as.vector(a*(y-h1t))
  m2<-cbind(cov_h0,z)*as.vector((1-a)*(h1t-h0c))
  m3<-cbind(cov_h0,z)*as.vector((1-a)*(y-h2c))
  m4<-cbind(cov_q0,w)*as.vector((1-a)*q0c-1)
  m5<-cbind(cov_q1,w)*(a*q0c*q1t-(1-a)*q0c)
  m6<-(a*q0c*q1t*(y-h1t)+(1-a)*q0c*(h1t-h0c)+h0c)-((1-a)*q0c*(y-h2c)+h2c)-psi[te_pos] 
  return(cbind(m1,m2,m3,m4,m5,m6))
}

p_mr_c<-function(data,h1_f,h0_f,q1_f,q0_f,w,z,trt,out){
  cov_h1<-model.matrix(h1_f,data=data);cov_h0<-model.matrix(h0_f,data=data)
  cov_q1<-model.matrix(q1_f,data=data);cov_q0<-model.matrix(q0_f,data=data)
  y=out;a=trt
  
  h1t_fit<-ivreg(y~-1+cov_h1+w|cov_h1+z,subset=a==1)
  h1t <- data.matrix(cbind(cov_h1,w))%*%h1t_fit$coefficients
  h0c_fit<-ivreg(h1t~-1+cov_h0+w|cov_h0+z,subset=a==0)
  h2c_fit<-ivreg(y~-1+cov_h0+w|cov_h0+z,subset=a==0)
  q0_data<-list(cov_q0=cov_q0,a=a,z=z,w=w);init<-rep(0,dim(cov_q0)[2]+dim(as.matrix(w))[2])
  q0c_fit<-gmm(moments_q0c, x = q0_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))
  q0c<-as.vector(q0_fun(gamma=coef(q0c_fit),int.cov_z=cbind(cov_q0,z),a=a))
  q1_data<-list(cov_q1=cov_q1,a=a,z=z,w=w,q0c=q0c);init<-rep(0,dim(cov_q1)[2]+dim(as.matrix(w))[2])
  q1t_fit<-gmm(moments_q1_t, x = q1_data, t0 = init, type = "iterative", crit = 1e-25, wmatrix = "optimal", method = "BFGS", control = list(reltol = 1e-25, maxit = 20000))

  init <- c(coef(h1t_fit),coef(h0c_fit),coef(h2c_fit),coef(q0c_fit),coef(q1t_fit),0);nparms=length(init)
  
  gmm_data<-list(y=y,a=a,z=z,w=w,cov_h1=cov_h1,cov_q1=cov_q1,cov_h0=cov_h0,cov_q0=cov_q0)
  
  gmm(p_mr_fun_c, x = gmm_data, t0 = init,
      eqConst=c(1:(nparms-1)),eqConstFullVcov=TRUE,
      crit = 1e-25,vcov="iid", wmatrix = "ident", optfct="nlminb",itermax=1000)
}

