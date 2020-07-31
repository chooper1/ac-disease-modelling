#fits parameters r_tilde, p, alpha, K_tilde for cumulative fatalities using Teismann's parametric model

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
  
  #parameters 
  r_tilde=par[1]
  p=par[2]
  alpha=par[3]
  K_tilde=par[4]
  
  
  #F is total fatalities
  dF_=r_tilde*F_^p*(1-(F_/K_tilde)^alpha)
  
  return(list(dF_))
}

#define the ODE system for fitting parameters to F and C as well as tau and mu_CFR
C_F_rates=function(t, x, par){
  with(as.list(c(x, par)), {
  #parameters
  r=par[1]
  p=par[2]
  alpha=par[3]
  K=par[4]
  mu_CFR=par[5]
  tau=par[6]
  
  r_tilde=r*(1-p)/(mu_CFR)
  K_tilde=mu_CFR*K
  
  
  #F is total fatalities
  dF_=r_tilde*F_^p*(1-(F_/K_tilde)^alpha)
  dC=r*C^p*(1-(C/K)^alpha)
  
  list(c(dF_, dC))
  })
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
  
  #parameter estimates returned from fit_param_F()
  F_parest=c(F_parest[1], F_parest[2], F_parest[3], F_parest[4])
  
  #separate values from before and after the first death (since ODE solver only starts at beginning of outbreak)
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
  
  #setting up the times to generate the curve after the outbreak starts
  times2=times2-start+1
  
  #we add value 1 at the beginning so the curve starts at the start of the deaths with the proper initial value
  times2=append(times2, 1, after=0)
  
  #generates the estimates at the required times
  F_out=ode(y=init, times=times2, func=F_rate, parms=F_parest)
  F_df=data.frame(F_out)
  times2=times2+start-1
  
  #relabels the times again so they are the same as the original
  colnames(F_df)=c("t", "pred_F")
  
  #removes the extra "1" we added to the times
  F_df=F_df[-1,]
  
  #builds the dataframe, adding the initial zeroes to the generated values
  y=c(y, F_df$pred_F)
  new_df=data.frame(times, y)
  
  return(new_df)
}

#input the times that have already been shifted by tau, and this function generates a fitted curve for the cumulative case data
generate_C=function(par, cases_C, times){

  #fits parameters to cumulative case data, starting at the first case
 start=min(which(cases_C>0, arr.ind=TRUE))
 cases_C=c(cases_C[start:length(cases_C)])
 fit_C=optim(par=par, fn=ssq2, cases=cases_C, control=list(parscale=c(1,1,1,10000)))
 C_parest=fit_C$par
 
 #removes times before the outbreak begins
 times2=times[-c(1:start-1)]
 init=cases_C[1]
 
 #solves the ODE for the given tiemes
 out=ode(y=init, times=times2, func=rate, parms=C_parest)
 outdf=data.frame(out)
 colnames(outdf)=c("times", "cases")
 
 #adds back in the "zero cases" at the beginning
 cases_C=rep(0, start-1)
 cases_C=append(cases_C, outdf$cases)
 df=data.frame(times, cases_C)
 
 return(df)
}

#for fitting tau and (optionally) mu_CFR
#par is a vector the guesses for the parameters to be fitted
#F_parest is the guesses for parameters for the fatality data, returned from fit_param_F()
ssq_C_F=function(par, cases_C, cases_F, F_parest){
  
  #sets default value for mu_CFR if only fitting for tau
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
  
  #generates F values for all the shifted times using the model parameters
  #these values have not been shifted, just taken from the curve given by the ODE solver
  F_df=generate_F(times_C, init, start, F_parest)
  
  #finds the absolute squared distance between the two curves
  ssqr=sum((F_df$y-C_df$cases_C*mu_CFR)^2)
  
  return(ssqr)
}

#this function is to visualize the effect of changing tau when fitting for tau and mu_CFR
#uses the automatic scaling of 0.05 by default for mu_CFR
#tau is a vector of tau values
plot_ssq_vs_tau=function(tau, cases_C, cases_F, F_parest){
  
  ssq=c()
  
  #calculates ssq for each tau value
  for(x in tau){
    ssq_tau=ssq_C_F(c(x), cases_C, cases_F, F_parest)
    ssq=append(ssq, ssq_tau)
  }
  
  df=data.frame(tau, ssq)
  
  #plots ssq vs. tau values
  plot=ggplot(data=df, aes(x=tau, y=ssq))+geom_line()
  print(plot)
}

#use the plot_shifted_scaled_cases() function from fit_datasets_F instead of this one
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


#for fitting tau, mu_CFR, and r, p, alpha, k all at once
#par is a vector the guesses for the parameters to be fitted
#F_parest is the guesses for parameters for the fatality data, returned from fit_param_F()
ssq_all_params_simultaneous=function(par, cases_C, cases_F){
  
  #parameters
  r=par[1]
  p=par[2]
  alpha=par[3]
  K=par[4]
  mu_CFR=par[5]
  tau=par[6]
  
  start_C=min(which(cases_C>0, arr.ind=TRUE))
  cases_C=c(cases_C[start_C:length(cases_C)])
  
  start_F=min(which(cases_F>0, arr.ind=TRUE))
  cases_F=c(cases_F[start_F:length(cases_F)])
  
  #times=seq(0, length(cases_C), 0.1)
  t=c(1:length(cases_C))
  #times=sort(union(times, t))
 # df=data.frame(t, cases_C)
  init=c(C=cases_C[1], F_=cases_F[1])
  
  #solves the ODE for times in t
  out_1=ode(y=init, times=t, func=C_F_rates, parms=par)
  
  cases_C_1=out_1[,"C"]
  cases_F_1=out_1[,"F_"]
  
  #shift times for cumulative cases by tau
  times_shifted=t+tau
  
  out_2=ode(y=init, times=times_shifted, func=C_F_rates, parms=par)
  
  cases_C_2=out_2[,"C"]
  
  #separate values from before and after the first death (since ODE solver only starts at beginning of outbreak)
  y=c()
  times2=c()
  for(x in times_shifted){
    if(x<start_F){
      y=append(y, 0)
    }
    else{
      times2=append(times2, x)
    }
  }

  #setting up the times to generate the fatality curve after the outbreak starts
  times2=times2-start_F+1
  
  #we add value 1 at the beginning so the curve starts at the start of the deaths with the proper initial value
  times2=append(times2, 1, after=0)
  
  #generates the estimates at the required times
  F_out=ode(y=init, times=times2, func=C_F_rates, parms=par)
  F_df=data.frame(F_out)
  
  #changes times_2 back to the correct values
  times2=times2+start_F-1
  
  #relabels the times again so they are the same as the original
  colnames(F_df)=c("t", "cases_F_2")
  
  #removes the extra "1" we added to the times
  F_df=F_df[-1,]
  
  #builds the dataframe, adding the initial zeroes to the generated values
  cases_F_2=c(y, F_df$cases_F_2)
  #new_df=data.frame(times_shifted, y)
  
  #finds the absolute squared distance between the death and cases data curves plus the distance for the two curves after the shift
  ssqr=sum(sqrt((cases_C_1-cases_C)^2)/mean(cases_C))+sum(sqrt((cases_F_1-cases_F)^2)/mean(cases_F))+sum((cases_F_2-cases_C_2*mu_CFR)^2)
  
  return(ssqr)
}
