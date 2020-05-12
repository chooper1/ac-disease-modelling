setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")
source("parameter_fitting.R")

#read data generated from SEIR SSA
SEIR_data<-read.csv("SSA runs tau=0.005.csv")


cases<-function(SEIR_data, run, tau){
  df=split(SEIR_data, f=SEIR_data$run)
  x=which(df[[run]][1]%%(1/tau)==1)
  totals=c()
  for(y in x){
    total=sum(df[[run]][y,5:14])
    totals=append(totals, total)
  }
  #df=data.frame(seq(1, length(totals)), totals)
  #colnames(df)=(c("day", "total"))
  #return(df)
  return(totals)
}

#fits r, p, alpha, and K for a given run
fit_param_SEIR=function(SEIR_data, run, tau){
  
  #format case data for a given region
  cases=cases(SEIR_data, run, tau)
  
  #starting guess for parameters
  par=c(r=2, p=1, alpha=1, K=5000)
  
  #performs the fit
  fit=optim(par=par, fn=ssq2, cases=cases, control=list(parscale=c(1,1,1,10000)))
  parest=fit$par
  
  return(parest)
}

fit_multiple_SEIR=function(SEIR_data, tau){
  df=data.frame(fit_param_SEIR(SEIR_data, 1, tau))
  colnames(df)="1"
  for(x in seq(2, length(split(SEIR_data, f=SEIR_data$run)))){
    print(x)
    params=data.frame(fit_param_SEIR(SEIR_data, x, tau))
    colnames(params)=x
    df=cbind(df, params)
    print(params)
  }
  return(df)
}

#write.csv(fitted, "4_param_fit_SEIR.csv", row.names=TRUE)
