import numpy as np
from distShape import distShape

# Simulates disease dynamics by a simple-minded individual-based algorithm,
# which is a state-, age- and time-after-infection dependent
# (discrete-time) stochastic process.

# syntax: [N,Cu,P,Ou,Y,x] = covid19_Net(bet,tau,ph,init,Net,T,sigma)
#
# input:  bet    (vector) infection probabilities:
#                  bet(1) = E
#                  bet(2) = I
#                  bet(3) = A
#         tau    (vector) disease-stage periods:
#                  tau(1) = mean incubation period (time from infection to
#                                                    max. of beta(x) )
#                  tau(2) = median time spent in E
#                  tau(3) = median time spent in I
#                  tau(4) = median time spent in J
#                  tau(5) = median time spent in A
#                  tau(6) = median time spent in H
#         ph     (vector) branching parameters
#                   ph(1) = symptomatic vs. asymptomatic
#                   ph(2) = mild vs. severe (hospital)
#                   ph(3) = recovery vs. death
#         init   initial number of infected (introduced into E+I_1 class)
#         Net    network (contact) matrix (should be quadratic, symmetric
#                and zero on the main diagonal; i.e., could be represented
#                by an upper-triangular matrix).
#                Net(j,jj) contains the average number of contacts between
#                individuals j and jj per day. A "contact" has to be defined
#                in some way; for instance,
#                  1 contact = "individuals j and jj are in the
#                               same room for 1 hour"
#                If the entry Net(j,jj) is not an integer, we interpret it as
#                  Net(j, jj) = "probability that j and jj are in the
#                                same room for 1h at any given day"
#                Since for a fixed day a contact either happens or it doesn't,
#                we flip an appropriately biased coin to make the decision.
#         T      number of days simulated
#         sigma  (vector)  shape paremters
#                  sigma(1)  beta
#                  sigma(2)  mu
#
# output: N      (Tx9 matrix) daily incidences
#                   N(:,1) =  newly infected
#                   N(:,2) =  new cases
#                   N(:,3) =  new asymptomatic cases
#                   N(:,4) =  new deaths
#                   N(:,5) =  newly recovered
#
#         Cu     (Tx5 matrix) cumulative incidences
#                   Cu(:,1) =  total infected
#                   Cu(:,2) =  total symptomatic cases
#                   Cu(:,3) =  total asymptomatic cases
#
#         P      (Tx5 matrix) daily prevalences
#                   P(:,1) =  E exposed (compartments E and I_1 of our compartment diagram)
#                   P(:,2) =  I symptomatic infected (compartment I_2)
#                   P(:,3) =  A asymptomatic infected
#                   P(:,4) =  J patients with mild symptoms
#                   P(:,5) =  H patients with severe symptoms in hospital
#                   P(:,6) =  S susceptibles
#
#         Ou     (Tx3 matrix) outcomes
#                   Ou(:,1) = R  total previously symptomatic recovered
#                   Ou(:,2) = F  total fatalities
#                   Ou(:,3) = B  total previously asymptomatic recovered
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

def covid19_Net(bet,tau,ph,init,Net,T,sigma,ndt):

    Npop = Net.shape[0];
    # adjust things to account for a time step that's different from "one day"
    dt = 1/ndt;   # time step
    Net = dt*Net; # since Net is given as contacts "per day", the number of contacts
                 # "per dt" is dt x C
    T = ndt*T;    # turn simulation time span T (in days) into number of iterations
                 # (index i)
    Y = np.zeros((Npop,1))
    x = (-1)*np.ones((Npop,1))
    #a = 90*np.random.uniform(0,1,size=(Npop,1))
    #Z = np.zeros((Npop,T))
    #y = np.zeros((Npop,T))

    N = np.zeros((T,5))
    Cu = np.zeros((T,3))
    P = np.zeros((T,6))
    Ou = np.zeros((T,3))

    # sprinkle the intial infectives over classes E and I_1
    drawind = np.random.randint(0, Npop, size=init)
    Y[drawind] = 1
    x[drawind] = np.random.randint(tau[1]) # update time-after-infection for initial

    # Fix R0!!!!!
    # to compute R0(t) on the fly:
    #cohortlength = np.maximum([np.ceil(tau[1])+np.ceil(tau[2]), np.ceil(tau[1])+np.ceil(tau[4])])
    #RR = np.zeros(Npop, Npop)  # matrix of infection events:
                               #   RR(j,jj) = "indiv. j has infected indiv. jj"
    R0 = np.zeros((T,1))
    #R0inf = np.zeros(Npop,1)   # vector keeping track of how many infections
                               # each individual causes
    #R0 =cell(1,2); #### Check this
    R0t = np.zeros((T,3))

    for k in range(0,T):

    #fix R0Theory
    #avecontacts = np.mean(np.sum(Net, axis=1))
    #xx = range(0,100,cohortlength)  #xx = linspace(0,cohortlength,100);
    #aveinfectiveness = (sum(beta(xx(xx<=tau(2)),1)) +ph(1)*sum(beta(xx(xx>tau(2) & xx<=tau(2)+tau(3)),2))+(1-ph(1))*sum(beta(xx(xx>tau(2) & xx<=tau(2)+tau(5)),3)))*100/cohortlength;
    #R0theory = aveinfectiveness*avecontacts;

        sig = sigma
        betas = beta(x,Y,sig,bet,tau)
        infind_ = betas > 10**(-3);
        b = betas[infind_]
        infind = np.nonzero(infind_)
        infind = infind[0]

        if len(b) > 0:
            infectives = len(infind) # number of infective infectives
            cNet = np.floor(Net[infind][:]) + np.floor(Net[infind][:]-np.floor(Net[infind][:])+np.random.uniform(0,1,size=(infectives,Npop)))
            cNet = cNet.astype(np.int64)
            B = np.repeat(b[:,np.newaxis], Npop, axis=1)
            YY = np.ones((infectives,1)) * np.transpose(Y)
            ZZ = np.zeros((infectives,Npop))
            yy = np.zeros((infectives,Npop))

            for n in range(1, np.max(cNet)+1):
                randseed = np.random.uniform(0,1,size=(infectives,Npop))
                randindex = np.logical_and(randseed < B, YY == 0)
                expr = np.logical_and(randindex, cNet >= n)
                ZZ[expr] = 1
                yy[expr] = 1

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
        

        # J -> R
        randseed = np.random.uniform(0,1,size=(len(jind),1))
        expr = randseed < mu(x[:,0][jind],4,sigma,tau,dt)
        expr_ind = np.nonzero(expr)[0]
        rind = jind[expr_ind]
        Y[rind] = 6

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

        N[k][4] = len(hrind) + len(rind)
        N[k][3] = len(dind)

        # A -> B
        randseed = np.random.uniform(0,1,size=(len(aind),1))
        expr = randseed < mu(x[:,0][aind],3,sigma,tau,dt)
        expr_ind = np.nonzero(expr)[0]
        bind = aind[expr_ind]
        Y[bind] = 8
                
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
             Cu[k][:] = Cu[k-1][:] + N[k][:3]
        else:
             Cu[k][:] = N[k][:3]

        #Fix R0 and R0t!!!

        # R0inf = sum(RR,2);
        # if i>cohortlength
        # R0num = sum(R0inf(x==cohortlength));
        # R0den = length(x(x==cohortlength));
        # if R0den >0
        #   R0time = R0num/R0den;
        #   R0t(i,1) = R0time;
        # end
        # R0t(i,2:3) = [R0num,R0den];
        # end
        #
        # R0{1,1} = R0theory;
        # R0{1,2} = RR;
        x[x > -1] += dt

    return [N,Cu,P,Ou,R0,R0t,Y,x]


def beta(x,n,sigma,bet,tau):
    b = np.zeros(len(x))
    if len(x) != len(n):
        return
    s = sigma[0]
    
    ind_1 = np.nonzero(n == 1)[0]
    ind_2 = np.nonzero(n == 2)[0]
    ind_3 = np.nonzero(n == 3)[0]

    b_1 = bet[0]*distShape(x,ind_1,tau[0],s,1)
    b_2 = bet[1]*distShape(x,ind_2,tau[0],s,1)
    b_3 = bet[2]*distShape(x,ind_3,tau[0],s,1)

    b = b_1 + b_2 + b_3;
    return b

def mu(x,n,sigma,tau,dt):
    s = sigma[1]
    taux3 = tau[1] + tau[2]
    taux4 = tau[1] + tau[3]
    taux5 = taux3 + tau[4]
    taux6 = taux3 + tau[5]
    if n == 1:    # E->I and E->A
       m = dt*distShape(x,range(0,len(x)),tau[1],s,2)
    elif n == 2:  # I->J and I->H
       m = dt*distShape(x,range(0,len(x)),taux3,s,2)
    elif n == 3:  # A->B
       m = dt*distShape(x,range(0,len(x)),taux4,s,2)
    elif n == 4:  # J->R
       m = dt*distShape(x,range(0,len(x)),taux5,s,2)
    elif n == 5:  # H->R and H->F
       m = dt*distShape(x,range(0,len(x)),taux6,s,2)
    else:
       m = dt*np.ones(len(x))
    return m

#def phi(x,n,ph):
#    if n == 1:      # symptomatic vs. asymptomatic
#       f = ph[0]
#    elif n == 2:    # mild vs. severe (hospital)
#       f = ph[1]
#    elif n == 3:    # recovery vs. death
#       f = ph[2]
#    else:
#       f = 0
#    return f
