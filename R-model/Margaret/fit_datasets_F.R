setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("parameter_fitting_F.R")



#read data from JHU
JHU_F_data <- read.csv("JHU_data/time_series_covid19_deaths_global.csv")
JHU_F_data<-t(JHU_F_data)
JHU_C_data<-read.csv("JHU_data/time_series_covid19_confirmed_global.csv")
JHU_C_data<-t(JHU_C_data)

#generates a vector of region labels for the estimates (JHU data)
regions=function(data){
  
  regions=c()
  for(x in seq(1, ncol(data))){
    #combines label from first and second row (region and country)
    regions=append(regions, paste(data[1, x], data[2, x]))
  }
  return(regions)
}

#fits r, p, alpha, and K for a given region
fit_param_F=function(region, data){

  #format case data for a given region
  cases=as.integer(data[5:nrow(data), region])
  cases=cases[!is.na(cases)]
  #discards data for days before initial outbreak
  start=min(which(cases>0, arr.ind=TRUE))
  cases=c(cases[start:length(cases)])
  
  #starting guess for parameters
  par=c(r_tilde=2, p=1, alpha=1, K_tilde=cases[length(cases)])
  
  #performs the fit
  fit=optim(par=par, fn=ssq_F, cases=cases, control=list(parscale=c(1,1,1,10)))
  parest=fit$par
  
  return(parest)
}

fit_multiple_F=function(data){
  
  #removes columns for countries with no cases
  few_cases=c()
  for(x in seq(1, ncol(data))){
    cases=as.integer(data[5:nrow(data), x])
    cases=cases[!is.na(cases)]
    if(cases[length(cases)]<=10){
      few_cases=append(few_cases, x)
    }
  }
  data=subset(data, select=-c(few_cases))
  
  #generates labels
  regions=regions(data)
  
  #performs the fit for the first element in the dataframe (assumption: the dataframe has at least 2 columns)
  paramdf=data.frame(fit_param_F(1, data))
  colnames(paramdf)=c(regions[1])
  
  #performs the fit for the rest of the regions, adds results to the dataframe
  #if all the data were working, the loop would be for seq(2, ncol(data))
  for(x in seq(2, ncol(data))){
    print(x)
    p=data.frame(fit_param_F(x, data))
    colnames(p)=c(regions[x])
    paramdf=cbind(paramdf, p)
    print(p)
  }
  return(paramdf)
}

fit_tau_mu_CFR=function(region, C_data, F_data){
  
  #estimates paramters for fatality data (discarding days before outbreak)
  F_parest=fit_param_F(region, F_data)
  
  #F_parest=c(F_parest[1], F_parest[2], F_parest[3], F_parest[4])
  
  #format case data for a given region
  cases_C=as.integer(C_data[5:nrow(C_data), region])
  cases_C=cases_C[!is.na(cases_C)]
  
  #format fatality data for a given region
  cases_F=as.integer(F_data[5:nrow(F_data), region])
  cases_F=cases_F[!is.na(cases_F)]

  #use this if just fitting tau
  #par=c(tau=20)
  # fit=optim(par=par, fn=ssq_C_F, cases_C=cases_C, cases_F=cases_F, F_parest=F_parest, method="Brent", lower=-30, upper=30)
  
  #for fitting tau and mu_CFR
  par=c(tau=20, mu_CFR=0.05)
  fit=optim(par=par, fn=ssq_C_F, cases_C=cases_C, cases_F=cases_F, F_parest=F_parest, control=list(parscale=c(1,.5)))
  
  parest=fit$par
  
  return(parest)
}
