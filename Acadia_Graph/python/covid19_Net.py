import numpy as np
from scipy.stats import norm
import math
from matrix_processor import matrix_processor


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
#         matrices (dictionary) dictionary of labelled contact matrices.
#                network (contact) matrices (should be quadratic, symmetric
#                and zero on the main diagonal; i.e., could be represented
#                by an upper-triangular matrix). This argument contains a list
#                of labelled matrices, with each corresponding to a specific
#                type of contact.
#                For each matrix, element [j,jj] contains the average number of
#                contacts between individuals j and jj per day. A "contact" has
#                to be defined in some way; for instance,
#                  1 contact = "individuals j and jj are in the
#                               same room for 1 hour"
#                If the entry M[j,jj] is not an integer, we interpret it as
#                  M[j,jj] = "probability that j and jj are in the
#                                same room for 1h at any given day"
#                Since for a fixed day a contact either happens or it doesn't,
#                we flip an appropriately biased coin to make the decision.
#         T      number of days simulated
#         sigma  (vector)  shape paremters
#                  sigma[0]  beta
#                  sigma[1]  mu
#         filters  (matrix of dictionaries)
#                filters for each day (to select and weight matrices)
#
# optional inputs:
#         quar      variable controlling quarantining. Default value: None
#                   if quar is of length 1 (eg. [s], where s is a scalar), then 
#                   an s-day quarantine is imposed. 
#                   if quar is of length 2 (eg. [s,t], where s,t are scalars), then 
#                   an s-day quarantine is imposed and a t-day lockdown is triggered
#                   on the first positive test result.
#
#         facility  variable controlling the testing facility. Default value: None
#                   a testing facility object can be passed in to trigger testing.
#
#         sampler   variable controlling the sampling strategy. Default value: None
#                   a sampling strategy object can be passed to trigger randomized
#                   testing strategies that aim to detect the virus prior to the 
#                   appearance of the first symptomatic case.
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
#                   Cu[:,4] =  total hospitalizations
#                   Cu[:,5] =  total tests returned
#                   Cu[:,6] =  total tests positive
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
#            xq1  (Npop x 1 vector) time in quarantine (index case)
#            xqc  (Npop x 1 vector) time in quarantine (secondary case)
#            xq   (Npop x 1 vector) time in quarantine (all cases)
#            xl   (scalar)          time in lockdown   

def covid19_Net(bet,tau,ph,init,matrices,T,sigma,filters,quar=None,facility=None,sampler=None):
    
    #initialize matrix processor object
    proc = matrix_processor()
    
    #initialize variables
    Npop = matrices['C'].shape[0]
    Y = np.zeros((Npop,1))
    x = (-1)*np.ones((Npop,1))
    N = np.zeros((T,9))
    Cu = np.zeros((T,7))
    P = np.zeros((T,6))
    Ou = np.zeros((T,3))

    #timer and temporary storage for implementing testing policies
    xq1 = (-1)*np.ones((Npop,1))
    xqc = (-1)*np.ones((Npop,1))
    xq = (-1)*np.ones((Npop,1))
    xl = -1

    # sprinkle the intial infectives over classes E and I_1
    drawind = np.random.randint(0, Npop, size=init)
    Y[drawind] = 1
    x[drawind] = np.random.randint(tau[1]) # update time-after-infection for initial

    cohortlength = np.maximum(math.ceil(tau[1])+math.ceil(tau[2]), math.ceil(tau[1])+math.ceil(tau[4]))
    #matrix of infection events:
    #RR[j,jj] = "indiv. j has infected indiv. jj"
    RR = np.zeros((Npop,Npop))
    R0inf = np.zeros((Npop,1))
    R0t = np.zeros((T,3))

    num_matrices = len(filters)
    rand_shift = np.random.randint(0, num_matrices)

    for k in range(0,T):
        NetInd = (k+rand_shift)%num_matrices
        if k == 0:
            quarind1 = []
            symptind = []
        #lockdown
        if quar is not None:
            locklength = quar[1]
            #campus-wide lockdown once first symptomatic case appears
            if (len(symptind) > 0 and xl == -1):
                Net = proc.process(Npop, matrices, filters[NetInd][1])
                xl = 0 # lockdown clock starts
            # reapply lockdown
            elif (xl > -1 and xl < locklength):
                Net = proc.process(Npop, matrices, filters[NetInd][1])
            # lockdown ends
            elif xl==locklength: # restore contacts with non-quarantined indiv.
                Net = proc.process(Npop, matrices, filters[NetInd][0])
                xl = -1
                # reset xl to allow for multiple lockdowns?
            else:
                Net = proc.process(Npop, matrices, filters[NetInd][0])
                
        else:
            Net = proc.process(Npop, matrices, filters[NetInd][0])

        #reapply quarantining restrictions to next day's matrix
        if quar is not None and facility is not None:
            quartime = quar[0] # length of quarantine period  (typ. 14 days)

            [releaseind,dum] = np.nonzero(xq1 == quartime)
            if len(releaseind) > 0:
                xq1[releaseind] = -1  #set clock back
                xq[releaseind] = -1   #set clock back

            # release secondary cases who are finished with quarantine
            [releaseind,dum] = np.nonzero(xqc == quartime)
            if len(releaseind) > 0:
                xqc[releaseind] = -1  #set clock back
                xq[releaseind] = -1   #set clock back

            #reapply quarantining
            reapplyind = np.array(xq > 0).flatten()
            Net[reapplyind] = 0
            Net[:,reapplyind] = 0

            # quarantine primary cases (who are not already quarantined)
            if len(quarind1) > 0:
                Net[quarind1] = 0
                Net[:,quarind1] = 0
                xq1[quarind1] = 0

            #so that classmates don't get simultaneously quarantined as a primary
            #and a secondary case
            xq = np.maximum(xq1, xqc)

            # quarantine all classmates of primary cases (who are not already
            # quarantined)
            Class = proc.process(Npop, matrices, {
                'C': {'weight': 1},
                })
            [dum1,dum2] = np.nonzero(Class[quarind1,:] > 0)
            quarindc = np.unique(dum2)
            if len(quarindc) > 0:
                quar_c = np.array(xq[quarindc] == -1)
                quar_c = quar_c.flatten()
                quarindc = quarindc[quar_c]
                Net[quarindc] = 0
                Net[:,quarindc] = 0
                xqc[quarindc] = 0
                if facility is not None:
                    #if we assume that all people in the same class as an index case are quarantined
                    Y_temp = Y[quarindc] > 0
                    Y_temp = np.array(Y_temp)
                    Y_temp = Y_temp.flatten()
                    facility.submit(k, quarindc, Y_temp)

            xq = np.maximum(xq1, xqc)  # quarantine clock for both primary and seconday cases

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
            x = np.add(x, np.reshape(yy_sum, (Npop,1)))
            N[k][0] = np.sum(x == 0)

        eind = np.nonzero(Y == 1)[0]
        iind = np.nonzero(Y == 2)[0]
        jind = np.nonzero(Y == 4)[0]
        hind = np.nonzero(Y == 5)[0]
        aind = np.nonzero(Y == 3)[0]

        #E -> I and E-> A
        randseed = np.random.uniform(0,1,size=(len(eind),2))
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][eind],1,sigma,tau), randseed[:,1] < ph[0])
        expr_ind = np.nonzero(expr)[0]
        symptind = eind[expr_ind]
        Y[symptind] = 2
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][eind],1,sigma,tau), randseed[:,1] >= ph[0])
        expr_ind = np.nonzero(expr)[0]
        asymptind = eind[expr_ind]
        Y[asymptind] = 3

        N[k][1] = len(symptind)
        N[k][2] = len(asymptind)

        # I -> J and I-> H
        randseed = np.random.uniform(0,1,size=(len(iind),2))
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][iind],2,sigma,tau), randseed[:,1] < ph[1])
        expr_ind = np.nonzero(expr)[0]
        mildind = iind[expr_ind]
        Y[mildind] = 4
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][iind],2,sigma,tau), randseed[:,1] >= ph[1])
        expr_ind = np.nonzero(expr)[0]
        hospind = iind[expr_ind]
        Y[hospind] = 3

        N[k][3] = len(mildind)
        N[k][4] = len(hospind)

        # J -> R
        randseed = np.random.uniform(0,1,size=(len(jind),1))
        expr = randseed < mu(x[:,0][jind],4,sigma,tau)
        expr_ind = np.nonzero(expr)[0]
        rind = jind[expr_ind]
        Y[rind] = 6

        N[k][5] = len(rind)

        # H -> R and H -> F
        randseed = np.random.uniform(0,1,size=(len(hind),2))

        expr = np.logical_and(randseed[:,0] < mu(x[:,0][hind],5,sigma,tau), randseed[:,1] < ph[2])
        expr_ind = np.nonzero(expr)[0]
        hrind = hind[expr_ind]
        Y[hrind] = 6
        expr = np.logical_and(randseed[:,0] < mu(x[:,0][hind],5,sigma,tau), randseed[:,1] >= ph[2])
        expr_ind = np.nonzero(expr)[0]
        dind = hind[expr_ind]
        Y[dind] = 7

        N[k][8] = len(hrind)
        N[k][6] = len(dind)

        # A -> B
        randseed = np.random.uniform(0,1,size=(len(aind),1))
        expr = randseed < mu(x[:,0][aind],3,sigma,tau)
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

        # Note:  we can add the use of the facility object to perform the tests
        # for folks in I.  When they go into I, some distribution(s) for compliance
        # to getting tested and probability of when they go.  When they get tested,
        # use facility.submit().  Then, we can pick them up in the getting of the
        # results a little lower.

        # Randomized testing with testing facilities
        if sampler is not None and facility is not None:
            # Determine who to "randomly" sample.  Skip dead and hospital
            subjects = sampler.sample(k, mask=(Y.isin([5,7])))
            # Consider excuding people who got tested today due to symptoms
            # Submit tests to the lab.  They'll be returned later
            if len(subjects) > 0:
                Y_temp = Y[subjects] > 0
                Y_temp = np.array(Y_temp)
                Y_temp = Y_temp.flatten()
                facility.submit(k, subjects, Y_temp)

        if facility is not None:
            # Get results that have come back from the lab today
            today_results_who, today_results = facility.results(k)

            # Okay, now what do we want to do with the positive test
            # results?  Add to a cumulative bin?
            if k == 1:
                Cu[k,5] = len(today_results_who)
                Cu[k,6] = np.sum(today_results)
            else:
                Cu[k,5] = Cu[k-1,5] + len(today_results_who)
                Cu[k,6] = Cu[k-1,6] + np.sum(today_results)

            if quar is not None:
                #if we assume that all people in the same class as an index case are quarantined
                quarind1 =  today_results_who[today_results == 1]
                neg_results = today_results_who[today_results == 0]
                sec_quar = np.nonzero(xqc > -1)
                quarind_release = [x for x in neg_results if x in sec_quar]
                #today_results_who[today_results == 0], np.nonzero(xqc > 0)
                #release all people in quarantine who test negative
                if len(quarind_release) > 0:
                    xqc[quarind_release] = -1
                    xq[quarind_release] = -1

        #increment all clocks
        x[x > -1] += 1
        xq1[xq1 > -1] += 1
        xqc[xqc > -1] += 1
        if xl > -1:
            xl += 1
        xq = np.maximum(xq1, xqc)

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

def mu(x,n,sigma,tau):
    s = sigma[1]
    taux3 = tau[1] + tau[2]
    taux4 = tau[1] + tau[3]
    taux5 = taux3 + tau[4]
    taux6 = taux3 + tau[5]
    if n == 1:    # E->I and E->A
       m = (1+np.tanh(s*(x-tau[1])))/2
    elif n == 2:  # I->J and I->H
       m = (1+np.tanh(s*(x-taux3)))/2
    elif n == 3:  # A->B
       m = (1+np.tanh(s*(x-taux4)))/2
    elif n == 4:  # J->R
       m = (1+np.tanh(s*(x-taux5)))/2
    elif n == 5:  # H->R and H->F
       m = (1+np.tanh(s*(x-taux6)))/2
    else:
       m = np.ones(len(x))
    return m
