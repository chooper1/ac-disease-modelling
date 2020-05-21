#fits parameters r, p, alpha, and K for Teismann's parametric model

#set to your working directory
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")

#load libraries
library(deSolve)
library(minpack.lm)
library(ggplot2)

#load data
#data <- read.csv("total_cases.csv")

#modify this line to choose region
#cases=data$NScases

#format data
#cases=cases[!is.na(cases)]
#t=c(1:length(cases))
#df=data.frame(t, cases)
#init=c(cases[1])


#define the ODE
rate=function(t, C, par){
 
  #parameters (alpha and K set)
  r=par[1]
  p=par[2]
  alpha=par[3]
  K=par[4]

  
  #c is total cases
  dC=r*C^p*(1-(C/K)^alpha)
  
  return(list(dC))
}


#function that calculates residuals (to be minimized) - generates a vector of residuals, use with nls.lm optimization
#ssqpar=function(par){
  
 # r=par[1]
  #p=par[2]
  #alpha=par[3]
  #K=par[4]

  #solves the ODE for times in t
#  out=ode(y=init, times=t, func=rate, parms=par)
  
  #formats predicted data from ODE
 # outdf=data.frame(out)
  #colnames(outdf)=c("t", "pred")
  
  #calculates residuals from ODE
#  ssqr=outdf$pred-df$cases
 # return(ssqr)
#}

#calculates the sum of squared residuals, use with optim optimization
ssq2=function(par, cases){
  
  par=c(r=par[1],p=par[2], alpha=par[3], K=par[4])
  
  times=seq(0, length(cases), 0.1)
  t=c(1:length(cases))
  times=sort(union(times, t))
  df=data.frame(t, cases)
  init=c(cases[1])
  
  #solves the ODE for times in t
  out=ode(y=init, times=times, func=rate, parms=par)
  
  #formats predicted data from ODE
  outdf=data.frame(out)
  colnames(outdf)=c("t", "pred")
  
  #calculates residuals from ODE
  ssqr=sum((outdf$pred[outdf$t %in% t]-df$cases)^2)
  return(ssqr)
}

#starting guess for parameters
par=c(r=2, p=1, alpha=1, K=5000)

#modified Levenberg-Marquardt algorithm to minimize residuals
#fitval=nls.lm(par=par, fn=ssqpar)

#minimize residuals with optim function
fit2=optim(par=par, fn=ssq2, cases=cases, control=list(parscale=c(1,1,1,10000)))


#plotting predicted and experimental

#simulating data based on estimated parameters
#parest=coef(fitval)
#parest2=fit2$par
#times=seq(min(t), max(t), 0.1)
#out=ode(y=init, times=times, func=rate, parms=parest2)
#outdf=data.frame(out)
#colnames(outdf)=c("t", "pred")

#formatting the plots
#plot=ggplot(data=outdf, aes(x=t, y=pred, color="red"))+geom_line()+geom_point(data=df, aes(x=t, cases, color="green"))+theme(legend.position="none")+labs(x="time (days)", y="Total cases")
#print(plot)
#plot=ggplot(data=outdf, aes(x=t, y=pred, color="red"))+geom_point(data=df, aes(x=t, cases, color="green"))+theme(legend.position="none")+labs(x="time (days)", y="Total cases")

#summary(fitval)
#parest2
