#set working directory
setwd("C:/Users/mjiho/ac-disease-modelling/")

#load libraries
library(ggplot2) #for plotting
library(reshape2) #for reshaping data
library(deSolve) #for solving DEs
library(minpack.lm) #least squares using lm algorithm

#load data
data_guangdon <- read.csv("R-model/Margaret/reformatted_data2.csv", header=TRUE, fileEncoding="UTF-8-BOM")

ggplot(data=data_guangdon, aes(x=Day, y=Confirmed))+geom_point(size=3)


rates=function(t, h, params){
 # define parameters
  c=params$c
  beta=params$beta
  q=params$q
  delta=params$delta
  lambda=params$lambda
  rho=params$rho
  delta_I=params$delta_I
  delta_q=params$delta_q
  gamma_I=params$gamma_I
  gamma_A=params$gamma_A
  gamma_H=params$gamma_H
  alpha=params$alpha
  theta=params$theta
  
  r=rep(0, length(h))
  r[1]=-(beta*c+c*q*(1-beta))*h["S"]*(h["I"]+theta*h["A"])+lambda*h["S_q"] 
  r[2]=beta*c*(1-q)*h["S"]*(h["I"]+theta*h["A"])-delta*h["E"]
  r[3]=delta*rho*h["E"]-(delta_I+alpha+gamma_I)*h["I"]
  r[4]=delta*(1-rho)*h["E"]-gamma_A*h["A"]
  r[5]=(1-beta)*c*q*h["S"]*(h["I"]+theta*h["A"])-lambda*h["S_q"]
  r[6]=beta*c*q*h["S"]*(h["I"]+theta*h["A"])-delta_q*h["E_q"]
  r[7]=delta_I*h["I"]+delta_q*h["E_q"]-(alpha+gamma_H)*h["H"] 
  r[8]=gamma_I*h["I"]+gamma_A*h["A"]+gamma_H*h["H"]
  
  return(list(r))
}

ICs=c(S=1000, E=0, I=200, A=0, S_q=0, E_q=0, H=0, R=0)

t=data_guangdon$Day

params=list(c=14.781, 
            beta=0.000000021, 
            q=0.000000189, 
            delta=0.142857143, 
            lambda=0.071428571, 
            rho=0.86834, 
            delta_I=0.13266, 
            delta_q=0.1259,
            gamma_I=0.33029,
            gamma_A=0.13978, 
            gamma_H=0.11624,
            alpha=0.0000178,
            theta=0.2)

out=ode(y=ICs, times=t, func=rates, parms=params)

ssq=function(params){
  #initial concentration
  ICs=c(S=1000, E=0, I=200, A=0, S_q=0, E_q=0, H=0, R=0)
  t=c(seq(0, 5, 62, 0.1), data_guangdon$Day)
  t=sort(unique(t))
  
  c=params[1]
  beta=params[2]
  q=params[3]
  delta=params[4]
  lambda=param[5]
  rho=params[6]
  delta_I=params[7]
  delta_q=params[8]
  gamma_I=params[9]
  gamma_A=params[10]
  gamma_H=params[11]
  alpha=params[12]
  theta=params[13]
  
  out=ode(y=ICs, times=t, func=rates, 
          parms=list(c, beta, q, delta, lambda, rho, delta_I, delta_q, gamma_I, gamma_A, gamma_H, alpha, theta))
  
  outdf=data.frame(out)
  outdf=outdf[outdf$Day %in% data_guangdon$Day,]
  
  preddf=melt(outdf, id.var="Day", value.name="Confirmed")
  expdf=melt(data_guangdon, id.var="Day", value.name="Confirmed")
  ssqres=preddf$Confirmed-expdf$Confirmed
  
  return(ssqres)
  
}