#code based on Dr. Watmough's code, model from Dr. Teismann, parameter values mostly from Teismann and Tang

library("ggplot2")


#set this to your working directory
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")

#source file
source("SEIAR_Teismann.R")

#read in parameter values
paramset<-read.csv("params_Teismann.csv")


#initial states
S0<-c(S=800, E=0, B=0, A=20, I_p=20, I_m=20, J=0, H=0, F=0, Q_m=0, Q_J=0, R=0)

#control reproduction number, when control measures are in place
#Rc=with(paramset, ((beta*rho*c*(1-q))/(delta_I+alpha+gamma_I)+(beta*c*theta*(1-rho)*(1-q))/gamma_A)*S0[1])

#Rc

#choose parameter sets
param=paramset[1, 2:15]


#final time for SSA
tf<-60

#output times for ODE solver
times<-seq(0, tf, by=0.01)

#number of simulations
A=repsim(S0, rates, nu, param=param, tf=360, simName, runs=20)


#plot E for every simulation
sampleruns<-ggplot(A, aes(x=t, y=E))+geom_line(aes(color=run))
print(sampleruns)

#plot L, I, A for one simulation
B=pluckrun(A, 3)
singlerun<-ggplot(B, aes(x=t, y=values))+geom_line(aes(color=ind))
print(singlerun)

res_ode<-ode(y=S0, times, func=rhs, param)

plot(res_ode[,1], res_ode[,'I_m'], col='red')

