setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("SEIR_fitting.R")



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

#fits r_tilde, p, alpha, and K_tilde for a given region
fit_param_F_SEIR=function(region, data){
  
  #format case data for a given region
  cases=as.integer(data[5:nrow(data), region])
  cases=cases[!is.na(cases)]
  #discards data for days before initial outbreak
  start=min(which(cases>0, arr.ind=TRUE))
  cases=c(cases[start:length(cases)])
  
  #starting guess for parameters
  par=c(r_tilde=2, p=1, alpha=1, K_tilde=cases[length(cases)])
  
  #performs the fit
  fit=optim(par=par, fn=ssq_F, cases=cases, control=list(parscale=c(1,1,1,10^floor(log10(cases[length(cases)])))))
  parest=fit$par
  
  return(parest)
}