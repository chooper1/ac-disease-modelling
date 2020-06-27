import numpy as np
import math
import scipy
from scipy.stats import norm

#Example for generating plots using this function:

#[R0theory,Norms,M,S,taud] = covid19_R0(bet,tau,ph,Net,T,sigma,R0t,N,x,RR)
#plt.plot(range(0,T), R0t[:,0], 'b')
#plt.plot(range(0,T), M, 'r')
#plt.plot(range(0,T), R0theory[0]*np.ones((T,1)), 'k')
#
#meancontacts = Norms[0]
#plt.title('R_0 = {}, ave. # contacts = {}'.format(R0theory[0], meancontacts))
#plt.legend(["R_t = N(t)/N(t-d)", "R_t from cohorts", "R_0 theory"])
#plt.show()


# syntax: [R0theory,Norms,M,S,taud] = covid19_R0(bet,tau,ph,Net,T,sigma,R0t,N,x,RR)
#
# input:  bet    (vector) infection probabilities:
#                  bet[0] = E
#                  bet[1] = I
#                  bet[2] = A
#
#         tau    (vector) disease-stage periods:
#                  tau[0] = mean incubation period (time from infection to
#                                                    max. of beta(x) )
#                  tau[1] = median time spent in E
#                  tau[2] = median time spent in I
#                  tau[3] = median time spent in J
#                  tau[4] = median time spent in A
#                  tau[5] = median time spent in H
#
#         ph     (vector) branching parameters
#                   ph[0] = symptomatic vs. asymptomatic
#                   ph[1] = mild vs. severe (hospital)
#                   ph[2] = recovery vs. death
#
#         Net    network (contact) matrices (should be quadratic, symmetric
#                and zero on the main diagonal; i.e., could be represented
#                by an upper-triangular matrix).
#                Net[i][j,jj] contains the average number of contacts between
#                individuals j and jj per day. A "contact" has to be defined
#                in some way; for instance,
#                  1 contact = "individuals j and jj are in the
#                               same room for 1 hour"
#                If the entry Net[i][j,jj] is not an integer, we interpret it as
#                  Net[i][j,jj] = "probability that j and jj are in the
#                                same room for 1h at any given day"
#                Since for a fixed day a contact either happens or it doesn't,
#                we flip an appropriately biased coin to make the decision.
#         T      number of days simulated
#         sigma  (vector)  shape paremters
#                  sigma[0]  beta
#                  sigma[1]  mu
#
#         R0T    output from covid19_Net.py
#         RR     output from covid19_Net.py
#         N      output from covid19_Net.py
#         x      output from covid19_Net.py
#
# output: R0theory (1 x 3 vector) 3 versions of R0
#                R0theory[0] = R0 computed using average number of contacts
#                R0theory[1] = R0 computed using dominant eigenvalue of Net
#                R0theory[2] = R0 computed using max. number of contacts
#         Norms  (1 x 3 vector) 3 characteristic numbers of Net
#                Norms[0] = average number of contacts
#                Norms[1] = dominant eigenvalue of Net
#                Norms[2] = max. number of contact occurring in the population
#         M  (T x 1 vector)  time series of "Rt from data"
#                M(i) = N(i)/N(i-d)
#         S      (Npop x Npop matrix)  matrix of "actual" serialintervals
#                  S[j,jj] = x[j] - x[jj]  if  RR[j,jj] = 1
#         taud   (scalar)  mean serial interval; taud = mean(S)

def covid19_R0(bet,tau,ph,Net,T,sigma,R0t,N,x,RR):
    Npop = Net.shape[0]
    dt = 1
    w = R0t.shape[0]
    t = np.linspace(1,T,w)
    M = np.zeros((T,1))

    S = np.zeros((Npop,Npop))
    # DE function to compute R0theory
    def Fde(Yde, x):
        Zde = np.zeros(3)
        Zde[0] = - mu(x,1,sigma,tau,dt)*Yde[0]
        Zde[1] = ph[0]*mu(x,1,sigma,tau,dt)*Yde[0] - mu(x,2,sigma,tau,dt)*Yde[1]
        Zde[2] = (1-ph[0])*mu(x,1,sigma,tau,dt)*Yde[0] - mu(x,3,sigma,tau,dt)*Yde[2]
        return Zde

    Xj = x * np.ones((1,Npop))
    Xjj = np.ones((Npop,1))*x.T
    S[RR == 1] = Xj[RR == 1] - Xjj[RR == 1]
    totalinf = np.sum(S>0)
    taud = np.sum(S) / totalinf
    ####################################
    #
    # Compute R0 according to the theory
    #
    ####################################
    cohortlength = np.maximum(math.ceil(tau[1])+math.ceil(tau[2]), math.ceil(tau[1])+math.ceil(tau[4]))
    maxcontacts = np.amax(Net)
    meancontacts = np.sum(Net) / Npop
    [eigs, eigv] = np.linalg.eig(Net)
    lambdaC = np.amax(eigs)

    xx = np.linspace(0, cohortlength, num=100)
    y0 = np.zeros(3)
    y0[0] = 1
    Rde = scipy.integrate.odeint(Fde, y0, xx)
    dummY = np.ones((len(xx)))
    R0th_ = np.sum(np.multiply(beta(xx,dummY,sigma,bet,tau),Rde[:,0])) + np.sum(np.multiply(beta(xx,2*dummY,sigma,bet,tau),Rde[:,1])) + np.sum(np.multiply(beta(xx,3*dummY,sigma,bet,tau),Rde[:,2]))
    R0th = R0th_ * cohortlength / 100

    R0theory = [R0th*meancontacts,R0th*lambdaC,R0th*maxcontacts]
    Norms = [meancontacts,lambdaC,maxcontacts]

    d = round(taud)
    d = d.astype(np.int64)
    for i in range(d, T-d):
        if (N[i-d,1]+N[i-d,2] > 0):
            M[i] = (N[i,1]+N[i,2])/(N[i-d,1]+N[i-d,2])

    return [R0theory,Norms,M,S,taud]

####################################
#
# function definitions
#
####################################
def beta(x,n,sigma,bet,tau):
    b = np.zeros(len(x))
    if len(x) != len(n):
        return
    s = sigma[0]

    ind_1 = np.nonzero(n == 1)[0]
    ind_2 = np.nonzero(n == 2)[0]
    ind_3 = np.nonzero(n == 3)[0]

    b_1 = np.zeros((len(x)))
    b_2 = np.zeros((len(x)))
    b_3 = np.zeros((len(x)))

    b_1[ind_1] = bet[0]*norm.pdf((x[ind_1] - tau[0]) / s) / s
    b_2[ind_2] = bet[1]*norm.pdf((x[ind_2] - tau[0]) / s) / s
    b_3[ind_3] = bet[2]*norm.pdf((x[ind_3] - tau[0]) / s) / s

    b = b_1 + b_2 + b_3;
    return b

def mu(x,n,sigma,tau,dt):
    s = sigma[1]
    taux3 = tau[1] + tau[2]
    taux4 = tau[1] + tau[3]
    taux5 = taux3 + tau[4]
    taux6 = taux3 + tau[5]
    if n == 1:    # E->I and E->A
       m = dt*(1+np.tanh(s*(x-tau[1])))/2
    elif n == 2:  # I->J and I->H
       m = dt*(1+np.tanh(s*(x-taux3)))/2
    elif n == 3:  # A->B
       m = dt*(1+np.tanh(s*(x-taux4)))/2
    elif n == 4:  # J->R
       m = dt*(1+np.tanh(s*(x-taux5)))/2
    elif n == 5:  # H->R and H->F
       m = dt*(1+np.tanh(s*(x-taux6)))/2
    else:
       m = dt*np.ones(len(x))
    return m
