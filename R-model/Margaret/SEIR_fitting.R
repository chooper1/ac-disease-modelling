library("deSolve")
library("zoo")
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("fit_datasets_F.R")


rhs_SEIR=function(t, y, par, fit, N, phi, beta_type){
  with(as.list(c(y, par, fit)), {
    
    if(t>=phi$times[1]){
    newdata=data.frame(times=as.numeric(t))
    phi=as.numeric(predict(fit, newdata=newdata))
    }
    else{
      phi=mean(phi$ratios)
    }
    
    if(beta_type=="equal"){
      beta_1=beta
      beta_2=beta
      beta_l=beta
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

ssq_SEIR=function(par, region, active_cases, cases_F, cases_R, mu_CFR=0.01, phi, times, start, fit, S0, pop, N, beta_type){
  
  if(beta_type=="different"){
    par=c(par[1], par[2], par[3], par[4], a=1,  
        kappa=1/1.2, eta=1/8, mu_E=1/4, mu_2=1/5)
  } else if(beta_type=="equal"){
    par=c(par[1], par[2], a=1, kappa=1/1.2, eta=1/8, mu_E=1/4, mu_2=1/5)
  }
  
    
  ode_soln=ode(y=S0, times, func=rhs_SEIR, par=par, fit=fit, N=N, phi=phi, beta_type=beta_type)
  
  I_2=ode_soln[,"I_2"]
  I_4=ode_soln[,"I_4"]
  R_2=ode_soln[,"R_2"]
  F_=ode_soln[,"F_"]
  
  ssq=sum(sqrt((active_cases-(I_2+I_4))^2))/mean(active_cases)+sum(sqrt((cases_F-F_)^2))/mean(cases_F)+sum(sqrt((cases_R-R_2)^2))/mean(cases_R)
  
  return(ssq)
  
}

fit_to_SEIR=function(region, C_data=JHU_C_data, F_data=JHU_F_data, R_data=JHU_R_data, mu_CFR=0.01, pop, beta_type="different", roll_size=1){
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
  cases_C=rollmean(cases_C, roll_size)
  times=rollmean(times, roll_size)
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  cases_F=c(cases_F[start:length(cases_F)])
  cases_F=rollmean(cases_F, roll_size)
  
  regions_R=regions(R_data)
  index_R=which(regions_R==region_name, arr.ind=TRUE)
  cases_R=as.integer(R_data[5:nrow(R_data), index_R])
  cases_R=cases_R[!is.na(cases_R)]
  cases_R=c(cases_R[start:length(cases_R)])
  cases_R=rollmean(cases_R, roll_size)
  
  active_cases=cases_C-cases_R
  
  S0=c(S=pop, E=0, L=0, I_1=0, I_2=active_cases[1], I_3=0, I_4=0, R_1=0, R_2=cases_R[1], F_=cases_F[1])
  N=pop
  
  if(beta_type=="different"){
  par=c(beta_l=0.1, beta_1=0.1, beta_2=0.1, rho=0.8)
  }
  else if(beta_type=="equal"){
    par=c(beta=0.1, rho=0.8)
  }
  
  #ODE_fit=optim(par=par, fn=ssq_SEIR, region=region, active_cases=active_cases, cases_F=cases_F, cases_R=cases_R, mu_CFR=mu_CFR, phi=phi, times=times, start=start, fit=fit, S0=S0, pop=pop, control=list(parscale=c(1,1, 1, 1, 1)))
  ODE_fit=optim(par=par, fn=ssq_SEIR, gr=NULL, region=region, active_cases=active_cases, cases_F=cases_F, cases_R=cases_R, mu_CFR=mu_CFR, phi=phi, times=times, start=start, fit=fit, S0=S0, pop=pop, N=N, beta_type=beta_type, method="L-BFGS-B", lower=c(0, 0, 0, 0, 0), upper=c(1, 1, 1, 100, 1))
  fit_par=ODE_fit$par
  
  return(fit_par)
  
}

#compartment is "C", "F", or "R"
#fit_param is output from fit_to_SEIR. If no parameters are passed, fit_to_SEIR will run
plot_SEIR_fit=function(region, C_data=JHU_C_data, F_data=JHU_F_data, R_data=JHU_R_data, mu_CFR=0.01, pop, beta_type="different",fit_par=NULL, compartment, roll_size=1){
  if(length(fit_par==0)){
    fit_par=fit_to_SEIR(region, C_data, F_data, R_data, mu_CFR, pop, beta_type, roll_size=roll_size)
  }
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
  cases_C=rollmean(cases_C, roll_size)
  times=rollmean(times, roll_size)
  C_df=data.frame(times, cases_C)
  
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]
  cases_F=c(cases_F[start:length(cases_F)])
  cases_F=rollmean(cases_F, roll_size)
  F_df=data.frame(times, cases_F)
  
  regions_R=regions(R_data)
  index_R=which(regions_R==region_name, arr.ind=TRUE)
  cases_R=as.integer(R_data[5:nrow(R_data), index_R])
  cases_R=cases_R[!is.na(cases_R)]
  cases_R=c(cases_R[start:length(cases_R)])
  cases_R=rollmean(cases_R, roll_size)
  R_df=data.frame(times, cases_R)
  
  active_cases=cases_C-cases_R
  
  S0=c(S=pop, E=0, L=0, I_1=0, I_2=active_cases[1], I_3=0, I_4=0, R_1=0, R_2=cases_R[1], F_=cases_F[1])
  N=pop
  
  
  if(beta_type=="different"){
    fit_par=c(fit_par[1], fit_par[2], fit_par[3], fit_par[4], a=1,  
          kappa=1/1.2, eta=1/8, mu_E=1/4, mu_2=1/5)
  } else if(beta_type=="equal"){
    fit_par=c(fit_par[1], fit_par[2], a=1, kappa=1/1.2, eta=1/8, mu_E=1/4, mu_2=1/5)
  }
  
  
  ode_soln=ode(y=S0, times, func=rhs_SEIR, par=fit_par, fit=fit, N=N, phi=phi, beta_type=beta_type)
  
  I_2=ode_soln[,"I_2"]
  I_4=ode_soln[,"I_4"]
  R_2=ode_soln[,"R_2"]
  F_=ode_soln[,"F_"]
  
  sim_cases_C=I_2+I_4
  sim_C_df=data.frame(times, sim_cases_C)
  sim_F_df=data.frame(times, F_)
  sim_R_df=data.frame(times, R_2)
  
  if(compartment=="C"){
    plot=ggplot(data=C_df, aes(x=times, y=cases_C))+geom_point()+geom_line(data=sim_C_df, aes(x=times, y=sim_cases_C))+labs(title=region_name)
    print(plot)
  }
  else if(compartment=="F"){
    plot=ggplot(data=F_df, aes(x=times, y=cases_F))+geom_point()+geom_line(data=sim_F_df, aes(x=times, y=F_))+labs(title=region_name)
    print(plot)
  }
  else if(compartment=="R"){
    plot=ggplot(data=R_df, aes(x=times, y=cases_R))+geom_point()+geom_line(data=sim_R_df, aes(x=times, y=R_2))+labs(title=region_name)
    print(plot)
  }
  
}
