function [EX,ES,EI,ER,P] = SIR_Mat(bet,tau,init,C,T,ndt,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Approximates probabilities for the SIR_Net model by a nonlinear matrix 
% iteration.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% syntax: [ES,EI,ER,P] = SIR_Mat(bet,tau,init,C,T,ndt,sigma)
%
% input:  bet    infection probability per contact
%
%         tau    (scalar or vector) 
%                   Markovian case:     (scalar) median time spent in I
%                   non-Markovian case: (2-vector)
%                       tau(1)   time from infection to maximal infectivity 
%                       tau(2)   median time spent in I
%                       tau(3)  (optional) tau-discretization 
%                               default tau(3) = 5 
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

markov = length(tau);  %toggles between Markovian and non-Markovian versions 
if markov == 3
    markov = 2; 
end     
%                  markov = 1      Markovian 
%                  markov = 2    non-Markovian 
if markov == 1
    tau = [tau,tau,1];
elseif length(tau) == 2
    tau = [tau 5];
end    



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
T = ndt*T;  % turn simulation time span T (in days) into number of iterations 
            % (index i)
            
Tau = T;  % 
ndtau = Tau*tau(3); 
dtau = 1/tau(3);             


%%%%%%%%%%%%%%%%%%%%%
% predefine variables
%%%%%%%%%%%%%%%%%%%%%

P = zeros(T,Npop,3);
P(1,:,1) = 1; 

if ndtau > 1
    Q = zeros(T,Npop,ndtau); 
end     

% sprinkle the intial infectives over class I
drawind = randi(Npop,init,1);   % draw individuals that will be infected 
                                % initially 
%drawprop = rand(init,1);  % draw initial probabilities for initially infected                         
%P(1,drawind,1) = drawprop;    % update the state
%P(1,drawind,2) = 1-drawprop;  % update state 

P(1,drawind,1) = 0;    
P(1,drawind,2) = 1; 

EX = zeros(T,Npop);


%
%avecontacts = mean(sum(C,2));
%xx = linspace(0,cohortlength,100);
%aveinfectiveness = sum(beta(xx(xx<=tau(2)),markov))*100/cohortlength;

%lambdaC = max(eig(C));
%lambdaC = 1;

%R0theory = [aveinfectiveness*avecontacts bet*lambdaC*tau(2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Actual algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:T-1

  %%%%%%%%%%%%%%%%%
  % define vector b
  %%%%%%%%%%%%%%%%%
 
  if markov == 1 
    v  = log(1-bet*P(i,:,2));  % auxiliary vector
    Mu = dt/tau(2);
  elseif markov == 2 
    x   =  linspace(1,Tau,ndtau);  % infection age   
    Bet = squeeze(P(i,:,2))*beta(x(1)) + (squeeze(Q(i,:,2:ndtau))*beta(x(2:ndtau))')'; 
    Mu  = dt*dtau*mu(x(1))*P(i,:,2); 
    v   = log(1-Bet);  % auxiliary vector 
  end 
  
  bm = exp(v*C);            % vector prod (1 - beta p ); i.e. b = 1 - bm 
  
  %%%%%%%%%%%%%%%%%%%%%%
  % update probabilities
  %%%%%%%%%%%%%%%%%%%%%%
 
  P(i+1,:,1) = bm.*P(i,:,1);
  P(i+1,:,2) = (1-Mu).*P(i,:,2)+P(i,:,1)-P(i+1,:,1); 
  P(i+1,:,3) = P(i,:,3)+Mu.*P(i,:,2);
  
  % extra probabilities to accommodate discretization of tau
  if ndtau >1
      Nutau  = dt*dtau*mu(x(2))*Q(i,:,1);  % the mu's associated with the Q's 
      Q(i+1,:,1) = Mu.*P(i,:,2)+(1-Nutau).*Q(i,:,1); 
      for ii = 2:ndtau
         Nutaum = Nutau; 
         Nutau = dt*dtau*mu(x(ii))*Q(i,:,ii);
         Q(i+1,:,ii) = Nutaum.*Q(i,:,ii-1)+(1-Nutau).*Q(i,:,ii); 
      end
      P(i+1,:,3) = P(i,:,3)+Nutau.*Q(i,:,ndtau);
  end
  
end 

%%%%%%%%%%%%%%%%%%%%%%%%%
% output expection values
%%%%%%%%%%%%%%%%%%%%%%%%%
     
ES = sum(P(:,:,1),2);
EI = sum(P(:,:,2),2);
  if ndtau > 1
      EI = EI + sum(Q,[2 3]);
  end    
ER = sum(P(:,:,3),2);


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