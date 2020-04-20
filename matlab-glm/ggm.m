function [Y] = ggm(par,t,init)
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

r = par(1);
p = par(2); 
T = max(t);
tn = length(t);

solODE = ode23(@ODEfun,[t(1),T],init); 
[C,N] = deval(solODE,t);
odeset('reltol',10^(-6)); 
Y(1:tn) = C;
Y(tn+1:2*tn) = N; 
   
   
   function z = ODEfun(t,y)
        % y = C
        z = r*y^p; 
   end   

    
end

