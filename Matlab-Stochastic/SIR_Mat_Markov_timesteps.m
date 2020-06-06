% set parameters 
init = 4; 
Npop = 1000;
C=300*rand(Npop,Npop)/Npop;
C = (C+C')/2;
C = C - diag(C).*eye(Npop); 

%https://www.cs.purdue.edu/homes/dgleich/demos/erdos_renyi/generate.html
%erdos-renyi graph
%n = Npop;
%g_prob = 0.25; %0.5 0.25
%rand('seed',100); % reseed so you get a similar picture
%G = rand(n,n) < g_prob;
%G = triu(G,1);
%G = G + G';
%C = G .* C;

%https://www.mathworks.com/help/matlab/math/build-watts-strogatz-small-world-graph-model.html
h = WattsStrogatz(Npop,7,0.1);
G = adjacency(h);
G = triu(G,1);
G = G + G';
C = G .* C;

lambdaC = max(eig(C)); 
T=100;
close all
kmark = 10;

m_values = [1,2,5,10];
beta = 0.1;
TAU = 10.1;

%%%%%%%%%%%%%%%%
%
% Markovian runs
%
%%%%%%%%%%%%%%%%
for m=1:4
    bet = beta*m_values(m);  
    for l=1:4
        tau = TAU * l/2;
        R0 = lambdaC*tau*bet; 
        figure(((m-1)*4+l));
        for j=1:4
            ndt = 3*j-2;
            t = linspace(1,T,ndt*T); 
            P = zeros(ndt*T,4);
            Ou = zeros(ndt*T,kmark);
            for k=1:kmark 
                [n,cu,p,ou,r0,r0t,y,x] = SIR_Net(bet,tau,init,C,T,ndt); 
                P(:,k) = p(:,1);
                Ou(:,k) = ou(:,1);
            end
            [ES,EI,ER] = SIR_Mat_Markov(bet,tau,init,C,T,ndt);
            
            %plotting exposed
            subplot(4,2,2*j-1) 
            hold;
            plot(t,P,'c');
            plot(t,EI,'b')
            hold off;
            
            %plotting cumulative cases
            subplot(4,2,2*j)
            hold;
            plot(t,Ou,'c');
            plot(t,ER,'b')
            hold off; 
        end   
        frame = getframe(figure(((m-1)*4+l)));
        im = frame2im(frame);
        s = strcat('SIR_Mat_Markov_Beta_', string(bet), '_Tau_', string(tau), '_WattsStrogatz.png');
        imwrite(im,s);
        (m-1)*4+l
    end
end    
