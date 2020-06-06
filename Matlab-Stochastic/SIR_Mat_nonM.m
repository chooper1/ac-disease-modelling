function [ES,EI,ER,P] = SIR_Mat_nonM(bet,tau,init,C,T,ndt,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Approximates probabilities for the SIR_Net model by a nonlinear matrix 
% iteration.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% syntax: [ES,EI,ER,P] = SIR_Mat_nonM(bet,tau,init,C,T,ndt,sigma)
%
% input:  bet    infection probability per contact
%
%         tau    (2 - vector)   characteristic time scales  
%                       tau(1)   time from infection to maximal infectivity 
%                       tau(2)   median time spent in I 
%
%         init   initial number of infected (introduced into I class)
%
%         C      contact (network) matrix (should be quadratic, symmetric
%                & zero on the main diagonal; i.e.; could be represented
%                by an upper-triangular matrix).
%                C(j,jj) contains the average number of contacts between
%                individuals j & jj per day. A "contact" has to be defined
%                in some way; for instance;
%                  1 contact = "individuals j & jj are in the
%                               same room for 1 hour"
%                If the entry C(j,jj) is not an integer, we interpret it as
%                  C(j, jj) = "probability that j & jj are in the
%                              same room for 1h at any given day"
%                Since for a fixed day a contact either happens or it doesn't,
%                we flip an appropriately biased coin to make the decision.
%
%         T      number of days simulated
%
%         sigma  (vector, optional)  shape paremters
%                  sigma(1)  beta(x) 
%                  sigma(2)  mu(x)
%                 default sigma =[1,1]
%
%         ndt    (integer, optional) number of time steps per day   
%
%
% output: ES      (T * ndt - vector) expectated values for S
%
%         EI      (T * ndt - vector) expectated values for I
%
%         ER      (T * ndt - vector) expectated values for R
%
%         P       (T*ndt x Npop x 3  - array) individual probabilities   
%                     P(i,j,k) = probability that individual j is in 
%                                compartment k at time t = i * dt 
%                 where 
%                         k = 1 corresponds to the S compartment
%                         k = 2 corresponds to the I compartment
%                         k = 3 corresponds to the R compartment
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Preliminary stuff
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = size(C);
Npop = v(1);


if nargin <= 5
    ndt = 1;
end    
if nargin <= 6 
    sigma = [1,1];
end 


% adjust things to account for a time step that's different from "one day"
dt = 1/ndt; % time step 
C = dt*C;   % since C is given as contacts "per day", the number of contacts 
            % "per dt" is dt x C

% infection-age vector and discretization 
x = linspace(1,T,ndt*T); 
ndtau = find(mu(x)>1-10^(-4),1); % number of I classes for discretizing 
                                 % beta(x) and mu(x); 10^(-4) is a threshold  

T = ndt*T;  % turn simulation time span T (in days) into number of iterations 
            % (index i)

   

%%%%%%%%%%%%%%%%%%%%%
% predefine variables
%%%%%%%%%%%%%%%%%%%%%

O = zeros(T,Npop);
O(1,:) = 1;  % initally everybody is susceptible 
P = zeros(T,Npop,ndtau);
R = zeros(T,Npop);
 
% sprinkle the intial infectives over class I
drawind = randi(Npop,init,1);   % draw individuals that will be infected 
                                % initially 
%drawprop = rand(init,1);  % draw initial probabilities for initially infected                         
%P(1,drawind,1) = drawprop;    % update the state
%P(1,drawind,2) = 1-drawprop;  % update state 

O(1,drawind) = 0;    
P(1,drawind,1) = 1; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Actual algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:T-1

  %%%%%%%%%%%%%%%%%
  % define vector b
  %%%%%%%%%%%%%%%%%
 
    Bet = (squeeze(P(i,:,:))*beta(x(1:ndtau))')'; 
    v   = log(1-Bet);  % auxiliary vector 
    bm = exp(v*C);     % vector prod (1 - beta p ); i.e. b = 1 - bm 
  
  %%%%%%%%%%%%%%%%%%%%%%
  % update probabilities
  %%%%%%%%%%%%%%%%%%%%%%
 
  O(i+1,:) = bm.*O(i,:);
  P(i+1,:,1) = O(i,:)-O(i+1,:); 
  for ii=2:ndtau 
      P(i+1,:,ii) = (1-dt*mu(x(ii-1)))*P(i,:,ii-1); 
  end 
  Mu  = (squeeze(P(i,:,:))*mu(x(1:ndtau))')'; 
  R(i+1,:) = R(i,:) + dt*Mu;
    
end 

%%%%%%%%%%%%%%%%%%%%%%%%%
% output expection values
%%%%%%%%%%%%%%%%%%%%%%%%%
     
ES = sum(O(:,:,1),2);
EI = sum(P(:,:,:),[2 3]);
ER = sum(R(:,:,1),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function definitions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % infection probability   
    function b = beta(x)   
          b = bet*distShape(x,tau(1),sigma(1),1);   
    end 

    % transition probabilties
    function m = mu(x)
          m = distShape(x,tau(2),sigma(2),2);  
    end   
        
end 