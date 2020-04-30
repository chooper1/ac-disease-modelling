#load libraries
library(deSolve)
library(minpack.lm)

#load and format data
NScases=c(21,28,41,51,68,73,90,110,122,127,147,170,193,207,236,262,293,310,342,373,407,428,445,474,517,547,579)
NSnewcases=c(7,13,10,17,5,17,20,12,5,20,23,23,14,29,26,31,17,32,31,34,21,17,29,43,30,32)
t=c(1:length(NScases))
df=data.frame(t, NScases)
init=c(NScases[1])

#define the ODE
rate=function(t, C, par){
 
  #parameters (alpha and K set)
  r=par[1]
  p=par[2]
  alpha=1
  K=10^8
  
  #c is total cases
  dC=r*(C^p)*(1-(C/K)^alpha)
  
  return(list(dC))
}


#function that calculates residuals (to be minimized)
ssqpar=function(par){
  
  r=par[1]
  p=par[2]
  alpha=1
  K=10^8
  
  #solves the ODE for times in t
  out=ode(y=init, times=t, func=rate, parms=par)
  
  #formats predicted data from ODE
  outdf=data.frame(out)
  colnames(outdf)=c("t", "pred")
  
  #calculates residuals from ODE
  ssqr=outdf$pred-df$NScases
  return(ssqr)
}

#starting guess for parameters
par=c(r=1.94, p=0.45)

#modified Levenberg-Marquardt algorithm to minimize residuals
fitval=nls.lm(par=par, fn=ssqpar)

summary(fitval)

