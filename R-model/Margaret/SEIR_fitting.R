library("deSolve")
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("fit_datasets_F.R")


#need to find a better way to find phi(t)
rhs_SEIR=function(t, y, par, phi){
  with(as.list(c(y, par)), {
  dS=-(1/(S+L+I_1+I_2+R_1+R_2))*S*(beta_l*L+beta_1*I_1+beta_2*I_2)
  dL=(1/(S+L+I_1+I_2+R_1+R_2))*S-a*kappa*L
  dI_1=a*(1-1/phi[which(phi$times==t), 2])*kappa*L-eta*I_1
  dI_2=a*(1/phi[which(phi$times==t), 2])*kappa*L-eta*I_2
  dR_1=eta*I_1
  dR_2=eta*I_2
  dF=eta*(1-rho)*I_2
  
  list(c(dS, dL, dI_1, dI_2, dR_1, dR_2, dF))
  })
}

ssq_SEIR=function(par, region, C_data, F_data, R_data, mu_CFR=1){

  
  par=c(beta_l=par[1], beta_1=par[2], beta_2=par[3], a=par[4], rho=par[5])
  
  phi=phi_vs_time(region, C_data, F_data, mu_CFR)
  tau_mu_CFR=fit_tau_mu_CFR(region, C_data, F_data)
  tau=tau_mu_CFR[1]
  
  phi$times=phi$times-tau
  times=phi$times
  start=times[1]
  
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  cases_C=c(cases_C[start:length(cases_C)])
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  cases_F=c(cases_F[start:length(cases_F)])
  
  cases_R=as.integer(R_data[5:nrow(R_data), region])
  cases_R=cases_R[!is.na(cases_R)]
  cases_R=c(cases_R[start:length(cases_R)])
  
  

  
  active_cases=cases_C-cases_R
  
  pop=971395
  s0=c(S=pop, L=0, I_1=(1-1/phi[1, 2])*active_cases[1], I_2=active_cases[1], R_1=(1-1/phi[1, 2])*cases_R[1], R_2=cases_R[1], F_=cases_F[1])
  
  ode_soln=ode(y=S0, times, func=rhs_SEIR, par=par, phi=phi)
  
  I_2=ode_soln[,"I_2"]
  R_2=ode_soln[,"R_2"]
  F_=ode_soln[,"F_"]
  
  ssq=sum((active_cases-R_2)^2)+sum((cases_F-F_)^2)+sum((cases_R-R_2)^2)
  
  return(ssq)
  
}

fit_to_SEIR=function(region, C_data, F_data, R_data, mu_CFR=1){
  par=c(beta_L=0.1, beta_1=0.1, beta_2=0.1, a=1, rho=0.8)
  
  fit=optim(par=par, fn=ssq_SEIR, region=region, C_data=C_data, F_data=F_data, R_data=R_data, mu_CFR=mu_CFR, control=list(parscale=c(1,1, 1, 1, 1)))
  

}