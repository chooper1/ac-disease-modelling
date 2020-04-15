library("GillespieSSA")

library("deSolve")

simName<-"SEI(A)R model Tang"

rates<-c("lambda*S_q", 
         "(1-beta)*c*q*S*(I+theta*A)", 
         "beta*c*(1-q)*S*(I+theta*A)", 
         "beta*c*q*S*(I+theta*A)", 
         "delta*(1-rho)*E", 
         "delta*rho*E", 
         "gamma_A*A", 
         "gamma_I*I", 
         "delta_I*I", 
         "delta_q*E_q", 
         "gamma_H*H", 
         "alpha*I", 
         "alpha*H")

nu<-matrix(c( 1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0,-1, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0,-1,-1, 0, 0,-1, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, 1,-1, 0,-1,
              0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0,
             -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0,-1, 0, 0, 0),
              nrow=8, byrow=TRUE)

rhs<-function(t, x, parms=NULL) {
  with(as.list(c(x, parms)),{
    
    dS=-(beta*c+c*q*(1-beta))*S*(I+theta*A)+lambda*S_q; 
    dE=beta*c*(1-q)*S*(I+theta*A)-delta*E;
    dI=delta*rho*E-(delta_I+alpha+gamma_I)*I; 
    dA=delta*(1-rho)*E-gamma_A*A;
    dS_q=(1-beta)*c*q*S*(I+theta*A)-lambda*S_q;
    dE_q=beta*c*q*S*(I+theta*A)-delta_q*E_q;
    dH=delta_I*I+delta_q*E_q-(alpha+gamma_H)*H; 
    dR=gamma_I*I+gamma_A*A+gamma_H*H;
    list(c(dS, dE, dI, dA, dS_q, dE_q, dH, dR))
  })
}

#control reproduction number

Rc=function(param){
  with(param, ((beta*rho*c*(1-q))/(delta_I+alpha+gamma_I)+(beta*c*theta*(1-rho)*(1-q))/gamma_A)*S0)
}

repsim=function(S0, rates, nu, param, tf, simName, runs, method=ssa.d()) {
  A=data.frame(NULL)
  for (run in seq(1, runs)) {
  res_ssa<-ssa(x0=S0,
                 a=rates, nu, parms=param, tf, 
                 method=method, simName, 
                 verbose=FALSE)
                 
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