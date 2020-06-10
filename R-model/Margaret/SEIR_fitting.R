library("deSolve")
library("zoo")
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("fit_datasets_F.R")


rhs_SEIR=function(t, y, par, fit, N, phi){
  with(as.list(c(y, par, fit)), {
    
    if(t>=phi$times[1]){
    newdata=data.frame(times=as.numeric(t))
    phi=as.numeric(predict(fit, newdata=newdata))
    }
    else{
      phi=mean(phi$ratios)
    }
    
    dS=-(S/N)*(beta_l*L+beta_1*I_1+beta_2*I_2)
    dE=(S/N)*(beta_l*L+beta_1*I_1+beta_2*I_2)-mu_E*E
    dL=mu_E*E-a*kappa*L
    dI_1=(1-1/phi)*kappa*L-mu_2*I_1
    #dI_1=0.5*a*kappa*L-eta*I_1
    dI_2=a*(1/phi)*kappa*L-mu_2*I_2
    #dI_2=0.5*a*kappa*L-eta*I_2
    dI_3=mu_2*I_1-eta*I_3
    dI_4=mu_2*I_2-eta*I_4
    dR_1=eta*I_3
    dR_2=eta*rho*I_4
    dF=eta*(1-rho)*I_4
    
    list(c(dS, dE, dL, dI_1, dI_2, dI_3, dI_4, dR_1, dR_2, dF))
  })
}

ssq_SEIR=function(par, region, active_cases, cases_F, cases_R, mu_CFR=0.01, phi, times, start, fit, S0, pop, N){
  
  par=c(par[1], par[2], par[3], par[4], a=1,  
        kappa=1/1.2, eta=1/8, mu_E=1/4, mu_2=1/5)
  
  ode_soln=ode(y=S0, times, func=rhs_SEIR, par=par, fit=fit, N=N, phi=phi)
  
  I_2=ode_soln[,"I_2"]
  I_4=ode_soln[,"I_4"]
  R_2=ode_soln[,"R_2"]
  F_=ode_soln[,"F_"]
  
  ssq=sum((active_cases-(I_2+I_4)^2)+sum((cases_F-F_)^2)+sum((cases_R-R_2)^2))
  
  return(ssq)
  
}

fit_to_SEIR=function(region, C_data=JHU_C_data, F_data=JHU_F_data, R_data=JHU_R_data, mu_CFR=0.01, pop){
  region_name=regions(C_data)[region]
  
  phi=phi_vs_time(region, C_data, F_data, mu_CFR)
  
  tau_mu_CFR=fit_tau_mu_CFR(region, C_data, F_data)
  tau=tau_mu_CFR[1]
  
  phi$times=phi$times-tau
  times=phi$times
  #start=times[1]
  
  fit=lm(ratios~poly(times, 4, raw=TRUE), data=phi)
  
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  start=min(which(cases_C>0, arr.ind=TRUE))
  times=c(start:length(cases_C))
  cases_C=c(cases_C[start:length(cases_C)])
  #cases_C=rollmean(cases_C, 5)
  #times=rollmean(times, 5)
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  cases_F=c(cases_F[start:length(cases_F)])
 # cases_F=rollmean(cases_F, 5)
  
  regions_R=regions(R_data)
  index_R=which(regions_R==region_name, arr.ind=TRUE)
  cases_R=as.integer(R_data[5:nrow(R_data), index_R])
  cases_R=cases_R[!is.na(cases_R)]
  cases_R=c(cases_R[start:length(cases_R)])
 # cases_R=rollmean(cases_R, 5)
  
  active_cases=cases_C-cases_R
  
  S0=c(S=pop, E=0, L=0, I_1=0, I_2=active_cases[1], I_3=0, I_4=0, R_1=0, R_2=cases_R[1], F_=cases_F[1])
  N=pop
  
  par=c(beta_l=0.1, beta_1=0.1, beta_2=0.1, rho=0.8)
  
  #ODE_fit=optim(par=par, fn=ssq_SEIR, region=region, active_cases=active_cases, cases_F=cases_F, cases_R=cases_R, mu_CFR=mu_CFR, phi=phi, times=times, start=start, fit=fit, S0=S0, pop=pop, control=list(parscale=c(1,1, 1, 1, 1)))
  ODE_fit=optim(par=par, fn=ssq_SEIR, gr=NULL, region=region, active_cases=active_cases, cases_F=cases_F, cases_R=cases_R, mu_CFR=mu_CFR, phi=phi, times=times, start=start, fit=fit, S0=S0, pop=pop, N=N, method="L-BFGS-B", lower=c(0, 0, 0, 0, 0), upper=c(1, 1, 1, 100, 1))
  
  return(ODE_fit)
  
}
