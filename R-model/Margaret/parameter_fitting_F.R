#fits parameters r, p, alpha, K for cumulative fatalities using Teismann's parametric model

#set to your working directory
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")

#load libraries
library(deSolve)
library(minpack.lm)
library(ggplot2)
library(dde)

#define the ODE
F_rate=function(t, F_, par){
  
  #parameters (alpha and K set)
  r_tilde=par[1]
  p=par[2]
  alpha=par[3]
  K_tilde=par[4]
  
  
  #F is total fatalities
  dF_=r_tilde*F_^p*(1-(F_/K_tilde)^alpha)
  
  return(list(dF_))
}

#calculates the sum of squared residuals, use with optim optimization
ssq_F=function(par, cases){
  
  par=c(r_tilde=par[1],p=par[2], alpha=par[3], K_tilde=par[4])
  
  times=seq(0, length(cases), 0.1)
  t=c(1:length(cases))
  times=sort(union(times, t))
  df=data.frame(t, cases)
  init=c(cases[1])
  
  #solves the ODE for times in t
  out=ode(y=init, times=times, func=F_rate, parms=par)
  
  #formats predicted data from ODE
  outdf=data.frame(out)
  colnames(outdf)=c("t", "pred")
  
  #calculates residuals from ODE
  ssqr=sum((outdf$pred[outdf$t %in% t]-df$cases)^2)
  return(ssqr)
}

ssq_C_F=function(par, cases_C, cases_F, F_parest){
  #mu_CFR=par[1]
  #tau=0
  
  tau=par[1]
  #mu_CFR=par[2]
  #default if just fitting tau
  mu_CFR=0.047
  r_tilde=F_parest[1]
  p=F_parest[2]
  alpha=F_parest[3]
  K_tilde=F_parest[4]
  
  times_C=c(1:length(cases_C))+tau
  C_df=data.frame(times_C, cases_C)
  
  start=min(which(cases_F>0, arr.ind=TRUE))
  
  init=c(cases_F[start])
  F_out=ode(y=init, times=times_C, func=F_rate, parms=F_parest)
  
  F_df=data.frame(F_out)
  colnames(F_df)=c("t", "pred_F")
  
  ssqr=sum((F_df$pred_F*mu_CFR-C_df$cases_C))
  
  return(ssqr)
}
