setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("parameter_fitting_F.R")



#read data from JHU
JHU_F_data <- read.csv("JHU_data/time_series_covid19_deaths_global.csv")
JHU_F_data<-t(JHU_F_data)

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
  par=c(r_tilde=2, p=1, alpha=1, K_tilde=5000)
  
  #performs the fit
  fit=optim(par=par, fn=ssq_F, cases=cases, control=list(parscale=c(1,1,1,10000)))
  parest=fit$par
  
  return(parest)
}

fit_multiple_F=function(data){
  
  #removes columns for countries with no cases
  no_cases=c()
  for(x in seq(1, ncol(data))){
    cases=as.integer(data[5:nrow(data), x])
    cases=cases[!is.na(cases)]
    if(sum(cases)==0){
      no_cases=append(no_cases, x)
    }
  }
  data=subset(data, select=-c(no_cases))
  
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