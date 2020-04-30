#code based on Dr. Watmough's code, model from Dr. Teismann, parameter values mostly from Teismann and Tang

library("ggplot2")


#set this to your working directory
setwd("C:/Users/mjiho/ac-disease-modelling/R-model/Margaret/")

#source file
source("SEIAR_Teismann.R")

#read in parameter values
paramset<-read.csv("params_Teismann.csv")


#initial states
S0<-c(S=1000, E=0, B=0, A=0, I_p=20, I_m=0, J=0, H=0, F=0, Q_m=0, Q_J=0, R=0)

#choose parameter sets
param=paramset[1, 2:17]

#control reproduction number, when control measures are in place
Rc=with(param, (mu_E/(mu_E+kappa))*((beta_p/(mu_1+kappa))+(phi*mu_1+beta_m)/((mu_1+kappa)*(mu_2+lambda))+((1-phi)*mu_1*beta_A)/((mu_1+kappa)*(mu_A+lambda))))

Rc

#final time for SSA
tf<-60

#output times for ODE solver
times<-seq(0, tf, by=0.01)

#number of simulations
A=repsim(S0, rates, nu, param=param, tf=360, simName, runs=20)


#plot E for every simulation
sampleruns<-ggplot(A, aes(x=t, y=R))+geom_line(aes(color=run))
print(sampleruns)

#plot L, I, A for one simulation
B=pluckrun(A, 3)
singlerun<-ggplot(B, aes(x=t, y=values))+geom_line(aes(color=ind))
print(singlerun)

res_ode<-ode(y=S0, times, func=rhs, param)

plot(res_ode[,1], res_ode[,'S'], col='red')

