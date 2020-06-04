library("deSolve")
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("fit_datasets_F.R")

#need to find a better way to find phi(t)
rhs_SEIR=function(t, y, par, fit){
  with(as.list(c(y, par, fit)), {
  
  newdata=data.frame(times=as.numeric(t))
  phi=as.numeric(predict(fit, newdata=newdata))

  dS=-(1/(S+L+I_1+I_2+R_1+R_2))*S*(beta_l*L+beta_1*I_1+beta_2*I_2)
  dL=(1/(S+L+I_1+I_2+R_1+R_2))*S*(beta_l*L+beta_1*I_1+beta_2*I_2)-a*kappa*L
  dI_1=a*(1-1/phi)*kappa*L-eta*I_1
  #dI_1=0.5*a*kappa*L-eta*I_1
  dI_2=a*(1/phi)*kappa*L-eta*I_2
  #dI_2=0.5*a*kappa*L-eta*I_2
  dR_1=eta*I_1
  dR_2=eta*rho*I_2
  dF=eta*(1-rho)*I_2
  
  list(c(dS, dL, dI_1, dI_2, dR_1, dR_2, dF))
  })
}

ssq_SEIR=function(par, region, cases_C, cases_F, cases_R, mu_CFR=0.01, phi, times, start, fit, S0, pop){
  
  par=c(par[1], par[2], par[3], par[4], par[5], 
        kappa=1/6, eta=1/12)
  
  ode_soln=ode(y=S0, times, func=rhs_SEIR, par=par, fit=fit)
  
  I_2=ode_soln[,"I_2"]
  R_2=ode_soln[,"R_2"]
  F_=ode_soln[,"F_"]
  
  ssq=sum((active_cases-R_2)^2)+sum((cases_F-F_)^2)+sum((cases_R-R_2)^2)
  
  return(ssq)
  
}

fit_to_SEIR=function(region, C_data=JHU_C_data, F_data=JHU_F_data, R_data=JHU_R_data, mu_CFR=0.01, pop){
  region_name=regions(C_data)[region]
  
  phi=phi_vs_time(region, C_data, F_data, mu_CFR)
  
  tau_mu_CFR=fit_tau_mu_CFR(region, C_data, F_data)
  tau=tau_mu_CFR[1]
  
  phi$times=phi$times-tau
  times=phi$times
  start=times[1]
  
  fit=lm(ratios~poly(times, 4, raw=TRUE), data=phi)
  
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  cases_C=c(cases_C[start:length(cases_C)])
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  cases_F=c(cases_F[start:length(cases_F)])
  
  regions_R=regions(R_data)
  index_R=which(regions_R==region_name, arr.ind=TRUE)
  cases_R=as.integer(R_data[5:nrow(R_data), index_R])
  cases_R=cases_R[!is.na(cases_R)]
  cases_R=c(cases_R[start:length(cases_R)])
  
  active_cases=cases_C-cases_R
  
  S0=c(S=pop, L=0, I_1=(1-1/phi[1, 2])*active_cases[1], I_2=active_cases[1], R_1=(1-1/phi[1, 2])*cases_R[1], R_2=cases_R[1], F_=cases_F[1])
  
  
  par=c(beta_l=0.1, beta_1=0.1, beta_2=0.1, a=1, rho=0.8)
  
  ODE_fit=optim(par=par, fn=ssq_SEIR, region=region, cases_C=cases_C, cases_F=cases_F, cases_R=cases_R, mu_CFR=mu_CFR, phi=phi, times=times, start=start, fit=fit, S0=S0, pop=pop, control=list(parscale=c(1,1, 1, 1, 1)))
  
  return(ODE_fit)

}
