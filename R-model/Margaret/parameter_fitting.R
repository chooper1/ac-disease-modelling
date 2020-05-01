#fits parameters r and p for Teismann's parametric model

#set to your working directory
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")

#load libraries
library(deSolve)
library(minpack.lm)
library(ggplot2)

#load data
data <- read.csv("total_cases.csv")

#modify this line to choose region
cases=data$NScases

#format data
cases=cases[!is.na(cases)]
t=c(1:length(cases))
df=data.frame(t, cases)
init=c(cases[1])

alpha=1
K=10^8

#define the ODE
rate=function(t, C, par){
 
  #parameters (alpha and K set)
  r=par[1]
  p=par[2]

  
  #c is total cases
  dC=r*(C^p)*(1-(C/K)^alpha)
  
  return(list(dC))
}


#function that calculates residuals (to be minimized)
ssqpar=function(par){
  
  r=par[1]
  p=par[2]

  
  #solves the ODE for times in t
  out=ode(y=init, times=t, func=rate, parms=par)
  
  #formats predicted data from ODE
  outdf=data.frame(out)
  colnames(outdf)=c("t", "pred")
  
  #calculates residuals from ODE
  ssqr=outdf$pred-df$cases
  return(ssqr)
}

#starting guess for parameters
par=c(r=1.94, p=0.45)

#modified Levenberg-Marquardt algorithm to minimize residuals
fitval=nls.lm(par=par, fn=ssqpar)
summary(fitval)

#plotting predicted and experimental

#simulating data based on estimated parameters
parest=coef(fitval)
times=seq(1, t[length(t)], 0.1)
out=ode(y=init, times=times, func=rate, parms=parest)
outdf=data.frame(out)
colnames(outdf)=c("t", "pred")

#formatting the plots
plot=ggplot(data=outdf, aes(x=t, y=pred, color="red"))+geom_line()+geom_point(data=df, aes(x=t, cases, color="green"))+theme(legend.position="none")+labs(x="time (days)", y="Total cases")
print(plot)
