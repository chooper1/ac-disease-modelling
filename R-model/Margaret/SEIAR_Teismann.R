library("GillespieSSA")

library("deSolve")

simName<-"SEIAR model Teismann"

rates<-c("c*(beta_p*I_p+beta_m*I_m+beta_A*A)*S/(S+E+B+A+I_p+I_m+J+H+F+Q_m+Q_J+R)",
         "gamma",
         "mu_E*E",
         "(1-phi)*mu_1*I_p",
         "mu_2*A",
         "phi*mu_1*I_p",
         "theta*mu_2*I_m",
         "(1-theta)*mu_2*I_m",
         "lambda*I_m",
         "lambda*J",
         "mu_2*Q_m",
         "mu_J*J",
         "mu_J*Q_J",
         "rho*mu_H*H",
         "(1-rho)*mu_H*H", 
         "kappa*E",
         "kappa*A",
         "kappa*I_p")

nu<-matrix(c(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             +1,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,
              0, 0, 0, 0,+1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0,+1, 0,+1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0,
              0,+1,-1,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,
              0, 0, 0, 0, 0,+1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0,+1, 0, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0,+1, 0, 0, 0, 0, 0,-1,-1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0,+1, 0,-1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1, 0,-1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,+1,+1,+1, 0, 0, 0, 0),
           nrow=12, byrow=TRUE)

rhs<-function(t, x, parms=NULL) {
  with(as.list(c(x, parms)),{
    
    dS=-c*(beta_p*I_p+beta_m*I_m+beta_A*A)*S/(S+E+B+A+I_p+I_m+J+H+F+Q_m+Q_J+R); 
    dE=c*(beta_p*I_p+beta_m*I_m+beta_A*A)*S/(S+E+B+A+I_p+I_m+J+H+F+Q_m+Q_J+R)+gamma-mu_E*E-kappa*E;
    dB=mu_2*A
    dA=gamma+(1-phi)*mu_1*I_p-mu_2*A-kappa*A;
    dI_p=gamma+mu_E*E-mu_1*I_p-kappa*I_p;
    dI_m=phi*mu_1*I_p-theta*mu_2*I_m-lambda*I_m;
    dJ=theta*mu_2*I_m-mu_J*J-lambda*J;
    dH=(1-theta)*mu_2*I_m-mu_H;
    dF=(1-rho)*mu_H*H;
    dQ_m=lambda*I_m-mu_2*Q_m;
    dQ_J=mu_2*Q_m+lambda*J-mu_J*Q_J;
    dR=mu_J*J+mu_J*Q_J+rho*mu_H*H;
    list(c(dS, dE, dB, dA, dI_p, dI_m, dJ, dH, dF, dQ_m, dQ_J, dR))
  })
}

#control reproduction number

#Rc=function(param){
 # with(param, ((beta*rho*c*(1-q))/(delta_I+alpha+gamma_I)+(beta*c*theta*(1-rho)*(1-q))/gamma_A)*S0)
#}

repsim=function(S0, rates, nu, param, tf, simName, runs, method=ssa.d()) {
  A=data.frame(NULL)
  for (run in seq(1, runs)) {
    res_ssa<-ssa(x0=S0,
                 a=rates, nu, parms=param, tf, 
                 method=method, simName, 
                 verbose=FALSE, ignoreNegativeState = TRUE)
    
    A=rbind(A, data.frame(res_ssa$data, run=as.factor(run)))
  }
  return(A)
}

pluckrun=function(data, sim, col=3:5) {
  wide=subset(data, run==sim)
  long=stack(wide[,col])
  long$t=rep(wide$t, times=length(col))
  return(long)
}