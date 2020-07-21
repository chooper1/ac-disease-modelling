import matplotlib.pyplot as plt
import numpy as np
from covid19_Net import covid19_Net
#from covid19_R0 import covid19_R0

import testing_facility as tf
import graph_ops
#from graph_ops import matrixCM

def main():
    # Define parameters
    bet= [0.008, 0.008, 0.008]
    tau = [5.3500,5.2000,5.0000,12.0000,12.0000,10.0000]
    ph = [0.5, 1, 1]
    T=100
    sigma = [1, 1]

    quar = [14,5]
    locklevel = 0.7
    facility = tf.testing_facility(0.03, 0.02, [0.0, 0.5, 0.45, 0.04, 0.01])

    print('start loading matrix')

    # contact matrix

    data = graph_ops.loadEdgeList("graph_data_8Dj1.csv")
    C,M,R,S = graph_ops.edgesAsMatrices(data)

    matrices = {
        'C' : C,
        'M' : M,
        'R' : R,
        'S' : S
        }

    print('done loading matrix')
    #contacts = 12/7*C + 1/420*M + 1/4*R + 1.0*S
    a = {
        'C': {'weight': 12/7},
        'M': {'weight': 1/420},
        'R': {'weight': 1/4},
        'S': {'weight': 1}
        }
    b = {
        'R': {'weight': locklevel*1/4},
        'S': {'weight': locklevel}
        }
    filter = {0:{0:a, 1:b}, 1:{0:a, 1:b}, 2:{0:a, 1:b}, 3:{0:a, 1:b}, 4:{0:a, 1:b}, 5:{0:b, 1:b}, 6:{0:b, 1:b}}
    #running the covid19_Net simulation
    Num = 50
    CuQ = np.zeros((5,Num))

    for qq in range(0,5):
        if qq == 0:
            quar = [14]
        elif qq == 1:
            quar = [12]
        elif qq == 2:
            quar = [10]
        elif qq == 3:
            quar = [8]
        elif qq == 4:
            quar = [6] 
        
        for run in range(0,Num):
            print(quar)
            print(run)
            [n,cu,p,ou,R0t,R0inf,RR,y,x] = covid19_Net(bet,tau,ph,matrices,T,sigma,filter,quar=quar,facility=facility)
            #print(cu)
            #print(cu[0])
            #print(cu[:,0])
            #print(np.max(cu[:,0]))
            CuQ[qq,run] = np.max(cu[:,0])
            #[N,Cu,P,Ou,R0,Y,X] = covid19_Net_5Q(bet,tau,ph,init,Netw,T,quar);
            #print(CuQ)
            
    #plotting
    CuQ_outbreaks = CuQ>10
    plot_x = [14, 12, 10, 8, 6]
    plot_y = np.sum(CuQ_outbreaks, axis = -1)
    plt.plot(plot_x, plot_y)
    plt.ylabel('Number of outbreaks > 10')
    plt.xlabel('Quarantine Length')
    plt.title('Number of Outbreaks')
    plt.show()

    #[R0theory,Norms,M,S,taud] = covid19_R0(bet,tau,ph,C[0],T,sigma,R0t,n,x,RR)

    #plt.plot(range(0,T), R0t[:,0], 'b')
    #plt.plot(range(0,T), M, 'r')
    #plt.plot(range(0,T), R0theory[0]*np.ones((T,1)), 'k')

    #meancontacts = Norms[0]
    #print(Norms)
    #plt.title('R_0 = {}, ave. # contacts = {}'.format(R0theory[0], meancontacts))

    #plt.legend(["R_t = N(t)/N(t-d)", "R_t from cohorts", "R_0 theory"])
    #plt.show()



if __name__== "__main__":
    main()
