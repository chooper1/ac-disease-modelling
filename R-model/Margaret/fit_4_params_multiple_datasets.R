setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("parameter_fitting.R")



#read data from JHU
JHU_data <- read.csv("time_series_19-covid-Confirmed_archived_0325.csv")
JHU_data<-t(JHU_data)

#generates a vector of region labels for the estimates (JHU data)
regions=function(tdata){
  
  regions=c()
  for(x in seq(1, ncol(tdata))){
    #combines label from first and second row (region and country)
    regions=append(regions, paste(tdata[1, x], tdata[2, x]))
  }
  return(regions)
}

#fits r, p, alpha, and K for a given region
fit_param=function(region){
  
  #format case data for a given region
  cases=as.integer(tdata[5:nrow(tdata), region])
  cases=cases[!is.na(cases)]
  #discards data for days before initial outbreak
  start=min(which(cases>0, arr.ind=TRUE))
  cases=c(cases[start:length(cases)])
  
  #starting guess for parameters
  par=c(r=2, p=1, alpha=1, K=5000)
  
  #performs the fit
  fit=optim(par=par, fn=ssq2, cases=cases, control=list(parscale=c(1,1,1,10000)))
  parest=fit$par
  
  return(parest)
}


fit_multiple=function(tdata){
  
  #removes columns for countries with no cases
  no_cases=c()
  for(x in seq(1, ncol(tdata))){
    cases=as.integer(tdata[5:nrow(tdata), x])
    cases=cases[!is.na(cases)]
    if(sum(cases)==0){
      no_cases=append(no_cases, x)
    }
  }
 # tdata=subset(tdata, select=-c(no_cases))
  
  #generates labels
  regions=regions(tdata)
  
  #performs the fit for the first element in the dataframe (assumption: the dataframe has at least 2 columns)
  paramdf=data.frame(fit_param(1))
  colnames(paramdf)=c(regions[1])
  
  #performs the fit for the rest of the regions, adds results to the dataframe
  #if all the data were working, the loop would be for seq(2, ncol(tdata))
  for(x in c(seq(2, 209), seq(401, 493))){
    print(x)
    p=data.frame(fit_param(x))
    colnames(p)=c(regions[x])
    paramdf=cbind(paramdf, p)
    print(p)
  }
  return(paramdf)
}

#write fitted parameters to csv, fit is the output from fit_multiple(tdata)
write.csv(fit, "4_param_fit.csv", row.names=TRUE)