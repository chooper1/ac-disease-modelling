function y = distShape(x,mu,sigma,s)
%y = 2-2*x*exp(-x^2)/sqr(pi)+ erf(x);  
if s == 1
  %x
  %mu
  y = exp(-(x-mu).^2/(2*sigma^2))/(sigma*sqrt(2*pi));
  %y
elseif s == 2  
  y = (1+tanh(sigma*(x-mu)))/2;
elseif s == 3  
 y = zeros(size(x));
 y(x>=mu-sigma/2 & x<=mu+sigma/2) = 1/sigma;
else 
  y = zeros(size(x));
end 
end

