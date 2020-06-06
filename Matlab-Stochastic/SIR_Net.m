function [N,Cu,P,Ou,R0,R0t,Y,x,avecont] = SIR_Net(bet,tau,init,C,T,ndt,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulates SIR disease dynamics by a simple-minded individual-based
% algorithm, which is a state-dependent (discrete-time) stochastic process.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% syntax: [N,Cu,P,Ou,R0,R0t,Y,x,avecont] = SIR_Net(bet,tau,init,C,T,ndt,sigma)
%
% input:  bet    infection probability per contact
%
%         tau    (scalar or vector) 
%                   Markovian case:     (scalar) median time spent in I
%                   non-Markovian case: (2-vector)
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
% output: N      (Tx2 matrix) daily incidences
%                   N(:,1) =  newly infected
%                   N(:,2) =  number of newly recovered
%         Cu     (Tx1 matrix) cumulative incidences
%                   Cu(:,1) =  total infected
%         P      (Tx2 matrix) daily prevalences
%                   P(:,2) =  S susceptibles
%                   P(:,1) =  I infected
%         Ou     (Tx1 matrix) outcomes
%                   Ou(:,1) = R  total recovered
%         Y      final state (see "variables")
%         x      final time-after-infection (infection age; see "variables")
%
% variables: Y   (Npop x 1 vector) state, where
%                         Y = 0 corresponds to the S compartment
%                         Y = 1 corresponds to the I compartment
%                         Y = 2 corresponds to the R compartment
%
%            t    (T-vector)        chronological time   (index i)
%            x    (Npop x 1 vector) infection age (time after infection) 
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
%                  markov = 1    Markovian 
%                  markov = 2    non-Markovian 
if markov == 1
    tau = [tau,tau];
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
            
%%%%%%%%%%%%%%%%%%%%%
% predefine variables
%%%%%%%%%%%%%%%%%%%%%
Y = zeros(Npop,1);
x = (-dt)*ones(Npop,1);

N = zeros(T,2);
Cu = zeros(T,1);
P = zeros(T,2);
Ou = zeros(T,1);

% sprinkle the intial infectives over classes I
drawind = randi(Npop,init,1);   % draw individuals that will be infected on
                                % day t=1
Y(drawind) = 1;                     % update the state
x(drawind) = dt*randi(floor(tau(2)),init,1);  % update time-after-infection 
                                           % for initially infected

% to compute R0[t] on the fly:
cohortlength = ceil(tau(2));
RR = zeros(Npop,Npop);     % matrix of infection events:
                           % RR(j,jj) = "indiv. j infected indiv. jj at
                           %             some point during the outbreak" 
                           % (see also vector RRinf below)       

R0 =cell(1,2);
R0t =zeros(T,3);


avecontacts = mean(sum(C,2));
xx = linspace(0,cohortlength,100);
aveinfectiveness = sum(beta(xx(xx<=tau(2)),markov))*100/cohortlength;

%lambdaC = max(eig(C));
lambdaC = 1;

R0theory = [aveinfectiveness*avecontacts bet*lambdaC*tau(2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Actual algorithm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:T

  %%%%%%%%%%%%%%%%%%%
  % infections happen
  %%%%%%%%%%%%%%%%%%%

 infind = find(beta(x,Y)>10^(-3));    % get indices of infectives 
 b = beta(x(infind),Y(infind));       % vector of infection  probabilities

 if ~isempty(b)
    infectives = length(infind);    % number of infectives
    cC = floor(C(infind,:)) + floor(C(infind,:)-floor(C(infind,:))+rand(infectives,Npop));
    avecont(i) = mean(sum(cC,2));
                                    % determine # of contacts infectives have
    B = b*ones(1,Npop);             % infectives x NpAop matrix with
                                    % columns of b vectors
    YY = ones(infectives,1)*Y';     % ... rows of Y' vectors
    ZZ = zeros(infectives,Npop);    % dummy matrix to keep track of new infections
    yy = zeros(infectives,Npop);    % ..... to update time-after-infection variable
    for n=1:max(cC,[],"all")
        randseed = rand(infectives,Npop);
        Ind = find(cC>=n);    % pick out the index pairs (infind[j],jj) for
                              % which contact matrix elements are >= n
        Newinfectionsind = Ind(randseed(Ind) < B(Ind) & YY(Ind) ==0);
                                        % flip (biased) coin to decide whether
                                        % infections occur
        ZZ(Newinfectionsind) = 1;  % register infections in dummy matrix
        yy(Newinfectionsind) = 1;  % start time-after-infection clock for
                                   % newly infected individuals

        RR(infind,:) = RR(infind,:) + ZZ; % update the matrix of who infected whom

    end
    Y = Y + sum(ZZ,1)';  % collapsing ZZ gives vector of new infections
    x = x + dt*sum(yy,1)';  % all slots corresponding to susceptibles are -dt
                            % so adding dt's will set them to zero

    N(i,1) = length(x(x==0));         % number of newly infected

  end

    %%%%%%%%%%%%%%%%%%%%
    % disease progresses
    %%%%%%%%%%%%%%%%%%%%

    iind = find(Y==1);

    % I -> R
    randseed = rand(length(iind),1);
    rind = iind(randseed<mu(x(iind),Y(iind)));
    Y(rind) = 2;
    N(i,2) = length(rind);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % output time series etc.
    %%%%%%%%%%%%%%%%%%%%%%%%%
    P(i,2) = length(Y(Y==0));
    P(i,1) = length(Y(Y==1));
    
    if i>1
       Cu(i,1) = Cu(i-1,1) + N(i,1);
       Ou(i,1) = Ou(i-1,1) + N(i,2);   
    else 
       Cu(i,1) = N(i,1);
       Ou(i,1) = N(i,2); 
    end 
    
    % vector of number of infections by individual
    RRinf = sum(RR,2);    % R0inf(j) = number of infections caused by 
                          % individual j over the course of the whole outbreak 

    % compute instantaneous (i.e. time-dependent) reproductive number                       
    if i>cohortlength
    R0num = sum(RRinf(x==cohortlength)); % number of infections caused by 
                                         % all individuals (i.e. by this 
                                         % "cohort" of infectives) whose 
                                         % infective period just ended 
    R0den = length(x(x==cohortlength));  % number of infectives(i.e. size 
                                         % of the cohort) whose infective
                                         % period just ended 
    if R0den >0
      R0time = R0num/R0den;  % average number of infections per individual  
                             % of said cohort 
      R0t(i,1) = R0time;    
    end
    R0t(i,2:3) = [R0num,R0den];
    end
    
    % fill cell array with R0-like stuff 
    R0{1,1} = R0theory;
    R0{1,2} = RR;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time-after-infection clock ticks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x(x>-dt) = x(x>-dt) + dt;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function definitions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % infection probability   
     function b = beta(x,y)   
      if markov == 1  
          b = zeros(size(x));
          b(y==1) = bet;
      elseif markov == 2 
          b = bet*distShape(x,tau(1),sigma(1),1);
      end    
    end 

    % transition probabilties
    function m = mu(x,y)
       if markov == 1  
          m = zeros(size(x));
          m(y==1) = dt/tau(2); % if dt = 1, mu = 1/tau(2) is the probability 
                               % of recovery on a given day; so the prob. 
                               % of recovery during a time interval of 
                               % length dt is dt x mu.  
       elseif markov == 2 
          m = dt*distShape(x,tau(2),sigma(2),2);
       end  
    end   
        
end 