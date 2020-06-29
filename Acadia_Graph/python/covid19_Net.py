import numpy as np
from scipy.stats import norm
import math

# Simulates disease dynamics by a simple-minded individual-based algorithm,
# which is a state-, age- and time-after-infection dependent
# (discrete-time) stochastic process.

# syntax: [N,Cu,P,Ou,R0, R0T, Y,x] = covid19_Net(bet,tau,ph,init,Net,T,sigma)
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
#         init   initial number of infected (introduced into E+I_1 class)
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
#         ndt    number of timesteps per day
#         C_ind  indices for each day (to select between matrices)
#
#         test   choice of testing policy
#                default: "0" - no testing policy
#                         "1" - simple quarantine of length q_len for symptomatic cases
#                         "2" - simple quarantine while people are in symptomatic bucket
#         q_len  quarantine length
#
#
# output: N      (Tx9 matrix) daily incidences
#                   N[:,0] =  newly infected
#                   N[:,1] =  new symptomatic cases
#                   N[:,2] =  new asymptomatic cases
#                   N[:,3] =  number of individuals new in compartment J
#                   N[:,4] =  number of new hospililizations
#                   N[:,5] =  number of patients with mild symptoms newly recovered
#                   N[:,6] =  number of new deaths
#                   N[:,7] =  number of individuals new in compartment B
#                   N[:,8] =  number of hospiltalized patients newly recovered
#
#         Cu     (Tx5 matrix) cumulative incidences
#                   Cu[:,0] =  total infected
#                   Cu[:,1] =  total symptomatic cases
#                   Cu[:,2] =  total asymptomatic cases
#                   Cu[:,3] =  total mild cases
#                   Cu[:,4] =  total hospililizations
#
#         P      (Tx6 matrix) daily prevalences
#                   P[:,0] =  E exposed (compartments E and I_1 of our compartment diagram)
#                   P[:,1] =  I symptomatic infected (compartment I_2)
#                   P[:,2] =  A asymptomatic infected
#                   P[:,3] =  J patients with mild symptoms
#                   P[:,4] =  H patients with severe symptoms in hospital
#                   P[:,5] =  S susceptibles
#
#         Ou     (Tx3 matrix) outcomes
#                   Ou[:,0] = R  total previously symptomatic recovered
#                   Ou[:,1] = F  total fatalities
#                   Ou[:,2] = B  total previously asymptomatic recovered
#
#         R0t  (T x 3 matrix)
#                               R0t[:,0] = Rt (effective R0 from simulation)
#                               R0t[:,1] = total number of infections caused
#                                          by cohort whose infectious period
#                                          ended on day t
#                               R0t[:,2] = size of that cohort
#         Rinf  (Npop x 1 vector) 
#                               Rinf[j] = number of infections caused by individual j
#         RR  (Npop x Npop matrix of zeros and ones)
#                               RR[j,jj] = 1  <=>  indiv. j infected indiv. jj
#
#         Y      final state (see "variables")
#         x      final time-after-infection
#
# variables: Y   (Npop x 1 vector) state, where
#                         Y = 0 corresponds to the S compartment
#                         Y = 1 corresponds to the E +I_1 compartment
#                         Y = 2 corresponds to the I_2 compartment
#                         Y = 3 corresponds to the A compartment
#                         Y = 4 corresponds to the J compartment
#                         Y = 5 corresponds to the H compartment
#                         Y = 6 corresponds to the R compartment
#                         Y = 7 corresponds to the F compartment
#                         Y = 8 corresponds to the B compartment
#            t    (T-vector)        chronological time   (index i)
#            a    (Npop x 1 vector) individual age       (index j)
#            x    (Npop x 1 vector) time after infection (index k)

def covid19_Net(bet,tau,ph,init,NetGrouped,T,sigma,ndt,C_ind,test=0,q_len=14):

    Npop = NetGrouped[0].shape[0];
    # adjust things to account for a time step that's different from "one day"
    dt = 1/ndt;   # time step
    NetGrouped = [dt*x for x in NetGrouped] # since Net is given as contacts "per day", the number of contacts
                 # "per dt" is dt x C

    T = ndt*T;    # turn simulation time span T (in days) into number of iterations
                 # (index i)
    Y = np.zeros((Npop,1))
    x = (-1)*np.ones((Npop,1))

    N = np.zeros((T,9))
    Cu = np.zeros((T,5))
    P = np.zeros((T,6))
    Ou = np.zeros((T,3))

    #timer and temporary storage for implementing testing policies
    timers = np.zeros((Npop,1))

    # sprinkle the intial infectives over classes E and I_1
    drawind = np.random.randint(0, Npop, size=init)
    Y[drawind] = 1
    x[drawind] = np.random.randint(tau[1]) # update time-after-infection for initial

    #cohortlength = double(max(ceil(tau(2))+ceil(tau(3)),ceil(tau(2))+ceil(tau(5))))
    cohortlength = np.maximum(math.ceil(tau[1])+math.ceil(tau[2]), math.ceil(tau[1])+math.ceil(tau[4]))
    #matrix of infection events:
    #RR[j,jj] = "indiv. j has infected indiv. jj"
    RR = np.zeros((Npop,Npop))
    R0inf = np.zeros((Npop,1))
    R0t = np.zeros((T,3))

    num_matrices = len(C_ind)
    rand_shift = np.random.randint(0, num_matrices)

    for k in range(0,T):
        NetInd = (k+rand_shift)%num_matrices
        Net = NetGrouped[C_ind[NetInd]]

        #testing policies
        if test == 1: #simple q_len-day quarantine
            timer = timers > 0
            timer = timer[:,0]
            if np.sum(timer) > 0:
                Net[timer] = 0
                Net[:, timer] = 0
        elif test == 2: #quarantine while in symptomatic bucket
            timer = timers > 0
            timer = timer[:,0]
            if np.sum(timer) > 0:
                Net[timer] = 0
                Net[:, timer] = 0

        sig = sigma
        infind_ = beta(x,Y,sig,bet,tau) > 10**(-6)
        
        b = beta(x[infind_].reshape((np.sum(infind_),1)),Y[infind_].reshape((np.sum(infind_),1)),sig,bet,tau)
        infind = np.nonzero(infind_)[0]

        if len(b) > 0:
            infectives = len(infind) # number of infective infectives
            cNet = Net[infind,:]
            B = b*np.ones((infectives,Npop))
            YY = np.ones((infectives,1)) * np.transpose(Y)
            ZZ = np.zeros((infectives,Npop))
            yy = np.zeros((infectives,Npop))
            
            randseed = np.random.uniform(0,1,size=(infectives,Npop))
            Ind_ = cNet > 0
            Ind = np.nonzero(Ind_)
            BB = np.multiply(np.log(1-B), cNet)
            BB = 1-np.exp(BB)
            expr = np.zeros((infectives,Npop), dtype=bool)
            expr[Ind] = expr[Ind] + np.logical_and(randseed[Ind] < BB[Ind], YY[Ind] + ZZ[Ind] == 0)#.reshape((infectives,Npop))

            Newinfectionsind = np.zeros((infectives,Npop), dtype=bool)
            Newinfectionsind[expr] = Ind_[expr]
            if np.sum(Newinfectionsind) != 0:
                ZZ[Newinfectionsind] = 1
                yy[Newinfectionsind] = 1

            # update the matrix of who infected whom
            RR[infind,:] = RR[infind,:] + ZZ
            ZZ_sum = np.sum(ZZ,0)
            yy_sum = np.sum(yy,0)
            Y = np.add(Y,np.reshape(ZZ_sum, (Npop,1)))
            x = np.add(x, np.reshape(dt*yy_sum, (Npop,1)))
            N[k][0] = np.sum(x == 0)

        eind = np.nonzero(Y == 1)[0]
        iind = np.nonzero(Y == 2)[0]
        jind = np.nonzero(Y == 4)[0]
        hind = np.nonzero(Y == 5)[0]
        aind = np.nonzero(Y == 3)[0]

        #E -> I and E-> A
        randseed = np.random.uniform(0,1,size=(len(eind),2))
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][eind],1,sigma,tau,dt), randseed[:,1] < ph[0])
        expr_ind = np.nonzero(expr)[0]
        symptind = eind[expr_ind]
        Y[symptind] = 2
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][eind],1,sigma,tau,dt), randseed[:,1] >= ph[0])
        expr_ind = np.nonzero(expr)[0]
        asymptind = eind[expr_ind]
        Y[asymptind] = 3

        N[k][1] = len(symptind)
        N[k][2] = len(asymptind)


        # I -> J and I-> H
        randseed = np.random.uniform(0,1,size=(len(iind),2))
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][iind],2,sigma,tau,dt), randseed[:,1] < ph[1])
        expr_ind = np.nonzero(expr)[0]
        mildind = iind[expr_ind]
        Y[mildind] = 4
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][iind],2,sigma,tau,dt), randseed[:,1] >= ph[1])
        expr_ind = np.nonzero(expr)[0]
        hospind = iind[expr_ind]
        Y[hospind] = 3

        N[k][3] = len(mildind)
        N[k][4] = len(hospind)

        # J -> R
        randseed = np.random.uniform(0,1,size=(len(jind),1))
        expr = randseed < mu(x[:,0][jind],4,sigma,tau,dt)
        expr_ind = np.nonzero(expr)[0]
        rind = jind[expr_ind]
        Y[rind] = 6

        N[k][5] = len(rind)

        # H -> R and H -> F
        randseed = np.random.uniform(0,1,size=(len(hind),2))

        expr = np.logical_and(randseed[:,0] < mu(x[:,0][hind],5,sigma,tau,dt), randseed[:,1] < ph[2])
        expr_ind = np.nonzero(expr)[0]
        hrind = hind[expr_ind]
        Y[hrind] = 6
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][hind],5,sigma,tau,dt), randseed[:,1] >= ph[2])
        expr_ind = np.nonzero(expr)[0]
        dind = hind[expr_ind]
        Y[dind] = 7

        N[k][8] = len(hrind)
        N[k][6] = len(dind)

        # A -> B
        randseed = np.random.uniform(0,1,size=(len(aind),1))
        expr = randseed < mu(x[:,0][aind],3,sigma,tau,dt)
        expr_ind = np.nonzero(expr)[0]
        bind = aind[expr_ind]
        Y[bind] = 8
        N[k][7] = len(bind)

        P[k][5] += np.count_nonzero(Y == 0)
        P[k][0] += np.count_nonzero(Y == 1)
        P[k][1] += np.count_nonzero(Y == 2)
        P[k][2] += np.count_nonzero(Y == 3)
        P[k][3] += np.count_nonzero(Y == 4)
        P[k][4] += np.count_nonzero(Y == 5)
        Ou[k][0] += np.count_nonzero(Y == 6)
        Ou[k][1] += np.count_nonzero(Y == 7)
        Ou[k][2] += np.count_nonzero(Y == 8)

        if k > 1:
             Cu[k][0:5] = Cu[k-1][0:5] + N[k][0:5]
        else:
             Cu[k][0:5] = N[k][0:5]

        R0inf = np.sum(RR,1)
        if k>cohortlength:
            cohort_ind = (x==cohortlength)
            R0num = np.sum(R0inf[cohort_ind[:,0]])
            R0den = len(x[cohort_ind[:,0]])
            if R0den > 0:
                R0time = R0num/R0den
                R0t[k,0] = R0time
            R0t[k, 1] = R0num
            R0t[k, 2] = R0den

        #testing policies
        if test == 1: #simple q_len-day quarantine
            #increment timers
            timer = timers > 0
            timer = timer[:,0]
            timers[timer] = timers[timer] + 1
            #print(timers[timer])
            #remove people from isolation when they leave
            timer = timers > q_len
            timer = timer[:,0]
            timers[timer] = 0
            #new quarantined
            timers[symptind] = 1
        elif test == 2: #quarantine while in symptomatic bucket
            #quarantine symptomatic cases
            timers[symptind] = 1
            #remove from quarantine when people leave the symptomatic bucket
            timers[mildind] = 0
            timers[hospind] = 0

        x[x > -1] += dt

    return [N,Cu,P,Ou,R0t,R0inf,RR,Y,x]


def beta(x,n,sigma,bet,tau):
    b = np.zeros(len(x))
    if len(x) != len(n):
        return
    s = sigma[0]

    ind_1 = np.nonzero(n == 1)[0]
    ind_2 = np.nonzero(n == 2)[0]
    ind_3 = np.nonzero(n == 3)[0]

    b_1 = np.zeros((len(x),1))
    b_2 = np.zeros((len(x),1))
    b_3 = np.zeros((len(x),1))
     
    if len(ind_1) != 0:
        b_1[ind_1] = bet[0]*norm.pdf((x[ind_1] - tau[0]) / s)/ s
    if len(ind_2) != 0:
        b_2[ind_2] = bet[1]*norm.pdf((x[ind_2] - tau[0]) / s)/ s
    if len(ind_3) != 0:
        b_3[ind_3] = bet[2]*norm.pdf((x[ind_3] - tau[0]) / s)/ s

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
