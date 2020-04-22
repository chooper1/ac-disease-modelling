function [Y] = glm_Prev_P(par,t,init,Cinit)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%
% Computes C(t) and N(t) according to (1a), (1b) if the project PDF document 
% with K =oo abd alpha =0 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% syntax: [Y] = ggm(r,p,t) 
%
% input:  par parameters (vector) where 
%            par(1) = r  (empirical) early growth rate
%            par(2) = p  scaling parameter (p < 1 polynomial, p = 1 exponential) 
%         t  time vector (days)
%         init initial condition C(0)
%
% output: Y(1:length(t)) = C            cummulative incidence (total cases)
%         Y(length+1:2*length(t)) = N   daily incidence (new cases)
%         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = max(t);
mu = 1/14; 
tn = length(t);
[g_temp, g_prev] = glm(par, t, Cinit);

solODE = ode23(@ODEfun,t,init); 

[C,N] = deval(solODE,t);
odeset('reltol',10^(-6)); 
Y(1:tn) = C;

   function z = ODEfun(t,y) 
       [C0,N0] = deval(g_prev,t);
       z = N0 - mu*y;
   end  
end
