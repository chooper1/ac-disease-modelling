library("deSolve")
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("SEIAR_Teismann.R")

#read in parameter values
paramset<-read.csv("params_Teismann.csv")

#choose parameter set
param=paramset[1, 2:17]

rhs_SEIR=function(t, y, par, parset){
  
  
  #dS=-c*(beta_p*I_p+beta_m*I_m+beta_A*A)*S/(S+E+B+A+I_p+I_m+J+H+F+Q_m+Q_J+R); 
  #dE=c*(beta_p*I_p+beta_m*I_m+beta_A*A)*S/(S+E+B+A+I_p+I_m+J+H+F+Q_m+Q_J+R)+gamma-mu_E*E-kappa*E;
  #dB=mu_A*A
  #dA=gamma+(1-phi)*mu_1*I_p-mu_A*A-kappa*A;
  #dI_p=gamma+mu_E*E-mu_1*I_p-kappa*I_p;
  #dI_m=phi*mu_1*I_p-theta*mu_2*I_m-lambda*I_m;
  #dJ=theta*mu_2*I_m-mu_J*J-lambda*J;
  #dH=(1-theta)*mu_2*I_m-mu_H;
  dF=(1-rho)*mu_H*H;
  #dQ_m=lambda*I_m-mu_2*Q_m;
  #dQ_J=mu_2*Q_m+lambda*J-mu_J*Q_J;
  #dR=mu_J*J+mu_J*Q_J+rho*mu_H*H;
  list(c(dS, dE, dB, dA, dI_p, dI_m, dJ, dH, dF, dQ_m, dQ_J, dR))
}

#not right yet, we need to use the hospitalization data
ssq_F=function(par, cases_F){
  
  par=c(rho=par[1])
  mu_H=0.076923
  
  times=seq(0, length(cases_F), 0.1)
  t=c(1:length(cases_F))
  times=sort(union(times, t))
  df=data.frame(t, cases_F)
  init=c(cases_F[1])
  
  #solves the ODE for times in t
  out=ode(y=init, times=times, func=rhs_SEIR, parms=par)
  
  #formats predicted data from ODE
  outdf=data.frame(out)
  colnames(outdf)=c("t", "pred")
  
  #calculates residuals from ODE
  ssqr=sum((outdf$pred[outdf$t %in% t]-df$cases_F)^2)
  return(ssqr)
}