#fits parameters r, p, alpha, K for cumulative fatalities using Teismann's parametric model

#set to your working directory
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("parameter_fitting.R")

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

#times are the times at which we want to generate our curve, so times-tau
generate_F=function(times, init, start, F_parest){
  F_parest=c(F_parest[1], F_parest[2], F_parest[3], F_parest[4])
  y=c()
  times2=c()
  for(x in times){
    if(x<start){
      y=append(y, 0)
    }
    else{
      times2=append(times2, x)
    }
  }
  times2=times2-start+1
  times2=append(times2, 1, after=0)
  F_out=ode(y=init, times=times2, func=F_rate, parms=F_parest)
  F_df=data.frame(F_out)
  times2=times2+start-1
  colnames(F_df)=c("t", "pred_F")
  F_df=F_df[-1,]
  y=c(y, F_df$pred_F)
  new_df=data.frame(times, y)
  return(new_df)
}

#input the times that have already been shifted by tau, and this function generates a fitted curve for the cumulative case data
generate_C=function(par, cases_C, times){

 start=min(which(cases_C>0, arr.ind=TRUE))
 cases_C=c(cases_C[start:length(cases_C)])
 fit_C=optim(par=par, fn=ssq2, cases=cases_C, control=list(parscale=c(1,1,1,10000)))
 C_parest=fit_C$par
 times2=times[-c(1:start-1)]
 init=cases_C[1]
 out=ode(y=init, times=times2, func=rate, parms=C_parest)
 outdf=data.frame(out)
 colnames(outdf)=c("times", "cases")
 cases_C=rep(0, start-1)
 cases_C=append(cases_C, outdf$cases)
 df=data.frame(times, cases_C)
 return(df)
}

ssq_C_F=function(par, cases_C, cases_F, F_parest){
  tau=par[1]
  if(length(par)==2){
    mu_CFR=par[2]
  }
  else{
    mu_CFR=0.05
  }
  
  r_tilde=F_parest[1]
  p=F_parest[2]
  alpha=F_parest[3]
  K_tilde=F_parest[4]
  
  #shift times for cumulative cases by tau
  times_C=c(1:length(cases_C))+tau
  C_df=data.frame(times_C, cases_C)
  
  start=min(which(cases_F>0, arr.ind=TRUE))
  init=cases_F[start]
  
  F_df=generate_F(times_C, init, start, F_parest)
  
  #print((F_df$y*mu_CFR-C_df$cases_C))
  ssqr=sum((F_df$y-C_df$cases_C*mu_CFR)^2)
  
  return(ssqr)
}

plot_ssq_vs_tau=function(tau, cases_C, cases_F, F_parest){
  ssq=c()
  for(x in tau){
    ssq_tau=ssq_C_F(c(x), cases_C, cases_F, F_parest)
    ssq=append(ssq, ssq_tau)
  }
  df=data.frame(tau, ssq)
  plot=ggplot(data=df, aes(x=tau, y=ssq))+geom_line()
  print(plot)
}


plot_shifted_cases=function(par, cases_C, cases_F, F_parest){
  tau=par[1]
  mu_CFR=par[2]
  
  r_tilde=F_parest[1]
  p=F_parest[2]
  alpha=F_parest[3]
  K_tilde=F_parest[4]
  
  #shift times for cumulative cases by tau
  times_C=c(1:length(cases_C))+tau
  C_df=data.frame(times_C, cases_C)
  
  start=min(which(cases_F>0, arr.ind=TRUE))
  init=cases_F[start]
  
  F_df=generate_F(times_C, init, start, F_parest)
  
  plot=ggplot(data=C_df, aes(x=times_C, y=cases_C*mu_CFR, color="red"))+geom_line()+geom_line(data=F_df, aes(x=times, y=y, color="green"))+geom_line()
  print(plot)
}
