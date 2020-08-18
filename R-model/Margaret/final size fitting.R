library("deSolve")
library("rootSolve")
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("fit_datasets_F.R")

rhs_SEIIAR=function(t, y, par, fit, N, phi){
  with(as.list(c(y, par, fit)), {
    
    #sets phi to the average for t before start, selects correct phi value if after start
    if(t>=phi$times[1]){
      newdata=data.frame(times=as.numeric(t))
      phi=as.numeric(predict(fit, newdata=newdata))
    }
    else{
      phi=mean(phi$ratios)
    }

    dS=-beta*(I_1+I_2+A)*S/N
    dE=beta*(I_1+I_2+A)*S/N-mu_E*E
    dI_1=mu_E*E-1.5*mu*I_1
    dI_2=1.5*(1/phi)*mu*I_1-mu*I_2
    dC=1.5*(1/phi)*mu*I_1
    dA=1.5*(1-(1/phi))*mu*I_1-mu*A
    dR_1=mu*A
    dF_=(1-rho)*mu
    dR_2=rho*mu
    
    list(c(dS, dE, dI_1, dI_2, dC, dA, dR_1, dF_, dR_2))
    
  })
}

ssq_SEIIAR=function(par, region, cases_C, cases_F, mu_IFR=0.01, phi, times, start, fit, pop, mu_CFR){
  
  p_pre_0=par[3]
  p_post_0=par[4]
  
  names=c("beta", "mu_E", "mu", "p_pre_0", "p_post_0", "rho")
  #par=c(par[1], par[2], par[3], par[4], par[5], rho=1-mu_CFR)
 # par=c(par[1], par[2], 0.2, par[3], par[4], rho=1-mu_CFR)
  par=c(par[1], 1/4, par[2], par[3], par[4], rho=1-mu_CFR)
  names(par)=names
  
  N=pop
  
  if(phi$times[1]>start){
    phi_start=mean(phi$ratios)
  }else{
    index=which(phi$times==start)
    phi_start=phi$ratios[index]
  }
  
  names=c("S", "E", "I_1", "I_2", "C", "A", "R_1", "F_", "R_2")
  S0=c(S=N, E=N*p_pre_0/2, I_1=N*p_pre_0/2, I_2=p_post_0*N/2, C=cases_C[1], A=p_post_0*N/2, R_1=0, F_=cases_F[1], R_2=0)
  names(S0)=names
  
  #p_post_0=(S0[5]+S0[4])/N

  ode_soln=ode(y=S0, times, func=rhs_SEIIAR, par=par, fit=fit, N=N, phi=phi)
  
  #I_2=ode_soln[,"I_2"]
  #R_2=ode_soln[,"R_2"]
  C=ode_soln[,"C"]
  F_=ode_soln[,"F_"]
  #total_cases=I_2+R_2+F_
  
  ssq=sum(sqrt((cases_C-C)^2))/mean(cases_C)+sum(sqrt((cases_F-F_)^2))/mean(cases_F)
  
  return(ssq)
  
}

fit_to_SEIIAR=function(region, C_data=JHU_C_data, F_data=JHU_F_data, mu_IFR=0.01, pop, intervention=1, final_size_guess=5000){
  
  region_name=regions(C_data)[region]
  
  phi=phi_vs_time(region, C_data, F_data, mu_IFR, final_size_guess = final_size_guess)
  
  tau_mu_CFR=fit_tau_mu_CFR(region, C_data, F_data)
  tau=tau_mu_CFR[1]
  mu_CFR=tau_mu_CFR[2]
  
  phi$times=phi$times-tau
  times=phi$times
  
  fit=lm(ratios~poly(times, 6, raw=TRUE), data=phi)
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  start=max(min(c(which(cases_C>0, arr.ind=TRUE))), intervention)
  times=c(start:length(cases_C))
  cases_C=c(cases_C[start:length(cases_C)])

  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  cases_F=c(cases_F[start:length(cases_F)])

  #par=c(beta=0.38, mu_E=1/4, mu=1/5, p_pre_0=0.0001, p_post=0.0001)
  #par=c(beta=0.38, mu_E=1/4, p_pre_0=0.0001, p_post_0=0.0001)
  par=c(beta=0.38, mu=1/5, p_pre_0=0.0001, p_post_0=0.0001)
  
  ode_fit=optim(par=par, fn=ssq_SEIIAR, gr=NULL, region=region, cases_C=cases_C, cases_F=cases_F, mu_IFR=mu_IFR, phi=phi, times=times, start=start, fit=fit, pop=pop, mu_CFR=mu_CFR, method="L-BFGS-B", lower=c(0, 0, 0), upper=c(1, 1, 1))
  
  fit_par=ode_fit$par
  
  return(fit_par)
  
}

plot_SEIIAR_fit=function(region, C_data=JHU_C_data, F_data=JHU_F_data, mu_IFR=0.01, pop, intervention=1, fit_par=NULL){
  
  
  
}

fn<-function(z){
  0.7795282*z+log(1-z)-log(1-0.01068354)
}
uniroot(fn, c(0, 1))