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
    a = 90*np.random.uniform(0,1,size=(Npop,1))
    Z = np.zeros((Npop,T))
    y = np.zeros((Npop,T))

    N = np.zeros((T,5))
    Cu = np.zeros((T,3))
    P = np.zeros((T,6))
    Ou = np.zeros((T,3))

    # sprinkle the intial infectives over classes E and I_1
    drawind = np.random.randint(0, Npop, size=init)
    for i in range(0,init):
        Y[drawind[i]] = 1 # update the state
        x[drawind[i]] = np.random.randint(tau[1]) # update time-after-infection for initial


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
        infind = []
        b = []
        for i in range(0,len(betas)):
            if betas[i] > 10**(-3):
                infind.append(i)
                b.append(betas[i])

        if len(b) > 0:
            infectives = len(infind) # number of infective infectives
            cNet = np.floor(Net[infind][:]) + np.floor(Net[infind][:]-np.floor(Net[infind][:])+np.random.uniform(0,1,size=(infectives,Npop)))
            cNet = cNet.astype(np.int64)
            B = np.ones((infectives,Npop))
            for i in range(0,len(b)):
                for j in range(0,Npop):
                    B[i][j] = b[i]
            YY = np.ones((infectives,1)) * np.transpose(Y)
            ZZ = np.zeros((infectives,Npop))
            yy = np.zeros((infectives,Npop))

            for n in range(1, np.max(cNet)+1):
                randseed = np.random.uniform(0,1,size=(infectives,Npop))
                Ind = []
                for i in range(0,len(cNet)):
                    for j in range(0,len(cNet[0])):
                        if cNet[i][j] >= n:
                            Ind.append((i,j))

                index_ = 0
                for (q,r) in Ind:
                    if randseed[q][r].item() < B[q][r].item() and YY[q][r].item() == 0:
                        ZZ[(q,r)] = 1
                        yy[(q,r)] = 1

            ZZ_sum = np.sum(ZZ,0)
            yy_sum = np.sum(yy,0)
            for i in range(0,len(Y)):
                Y[i] = Y[i] + ZZ_sum[i]  # collapsing ZZ gives vector of new infections
                x[i] = x[i] + dt*yy_sum[i]

            N[k][0] = 0
            for i in range(0,len(x)):
                if x[i] == 0:
                    N[k][0] = N[k][0] + 1 # number of newly infected

        eind = []
        iind = []
        jind = []
        hind = []
        aind = []

        for i in range(0,len(Y)):
            if Y[i] == 1:
                eind.append(i)
            elif Y[i] == 2:
                iind.append(i)
            elif Y[i] == 4:
                jind.append(i)
            elif Y[i] == 5:
                hind.append(i)
            elif Y[i] == 3:
                aind.append(i)

        #E -> I and E-> A
        randseed = np.random.uniform(0,1,size=(len(eind),2))
        symptind = []
        asymptind = []
        for i in range(0,len(eind)):
            if randseed[i][0] < mu(x[eind[i]],1,sigma,tau,dt) and randseed[i][1] < ph[0]: #phi(a[eind[i]],1,ph):
                symptind.append(eind[i])
                Y[eind[i]] = 2
            elif randseed[i][0] < mu(x[eind[i]],1,sigma,tau,dt) and randseed[i][1] >= ph[0]: #phi(a[eind[i]],1,ph):
                asymptind.append(eind[i])
                Y[eind[i]] = 3

        N[k][1] = len(symptind)
        N[k][2] = len(asymptind)

        # I -> J and I-> H
        randseed = np.random.uniform(0,1,size=(len(iind),2))
        mildind = []
        hospind = []
        for i in range(0,len(iind)):
            if randseed[i][0] < mu(x[iind[i]],2,sigma,tau,dt) and randseed[i][1] < ph[1]: #phi(a[iind[i]],2,ph):
                mildind.append(iind[i])
                Y[iind[i]] = 4
            elif randseed[i][0] < mu(x[iind[i]],2,sigma,tau,dt) and randseed[i][1] >= ph[1]: #phi(a[iind[i]],2,ph):
                hospind.append(iind[i])
                Y[iind[i]] = 5

        # J -> R
        randseed = np.random.uniform(0,1,size=(len(jind),1))
        rind = []
        for i in range(0,len(jind)):
            if randseed[i] < mu(x[jind[i]],4,sigma,tau,dt):
                rind.append(jind[i])
                Y[jind[i]] = 6

        # H -> R and H -> F
        randseed = np.random.uniform(0,1,size=(len(hind),2))

        hrind = []
        dind = []
        for i in range(0,len(hind)):
            if randseed[i][0] < mu(x[hind[i]],5,sigma,tau,dt) and randseed[i][1] < ph[2]: #phi(a[hind[i]],3,ph):
                hrind.append(hind[i])
                Y[hind[i]] = 6
            elif randseed[i][0] < mu(x[hind[i]],5,sigma,tau,dt) and randseed[i][1] >= ph[2]: #phi(a[hind[i]],3,ph):
                dind.append(hind[i])
                Y[hind[i]] = 7

        N[k][4] = len(hrind) + len(rind)
        N[k][3] = len(dind)

        # A -> B
        randseed = np.random.uniform(0,1,size=(len(aind),1))
        bind = []
        for i in range(0,len(aind)):
            if randseed[i] < mu(x[aind[i]],3,sigma,tau,dt):
                bind.append(aind[i])
                Y[aind[i]] = 8

        for i in range(0,len(Y)):
            if Y[i] == 0:
                P[k][5] = P[k][5] + 1
            elif Y[i] == 1:
                P[k][0] = P[k][0] + 1
            elif Y[i] == 2:
                P[k][1] = P[k][1] + 1
            elif Y[i] == 3:
                P[k][2] = P[k][2] + 1
            elif Y[i] == 4:
                P[k][3] = P[k][3] + 1
            elif Y[i] == 5:
                P[k][4] = P[k][4] + 1
            elif Y[i] == 6:
                Ou[k][0] = Ou[k][0] + 1
            elif Y[i] == 7:
                Ou[k][1] = Ou[k][1] + 1
            elif Y[i] == 8:
                Ou[k][2] = Ou[k][2] + 1

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

        for i in range(0,len(x)):
            if x[i] > -1:
                x[i] = x[i] + dt

    return [N,Cu,P,Ou,R0,R0t,Y,x]


def beta(x,n,sigma,bet,tau):
    b = np.zeros(len(x))
    if len(x) != len(n):
        return
    s = sigma[0]
    ind_1 = []
    ind_2 = []
    ind_3 = []
    for i in range(0,len(x)):
        if n[i] == 1:
            ind_1.append(i)
        elif n[i] == 2:
            ind_2.append(i)
        elif n[i] == 3:
            ind_3.append(i)

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
