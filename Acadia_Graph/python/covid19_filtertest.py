import matplotlib.pyplot as plt
import numpy as np
from covid19_Net import covid19_Net
from covid19_R0 import covid19_R0

import testing_facility as tf
import graph_ops
#from graph_ops import matrixCM

def main():
    # Define parameters
    bet= [0.01, 0.01, 0.01]
    tau = [5, 6, 4, 8, 17, 10]
    ph = [0.5, 1, 1]
    init = 1
    T=100
    Npop = 4000
    sigma = [1, 1]

    quar = [14,3]
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
    [n,cu,p,ou,R0t,R0inf,RR,y,x] = covid19_Net(bet,tau,ph,init,matrices,T,sigma,filter,quar = quar,facility=facility)

    #plotting
    plot_p = np.zeros(T)
    plot_cu = np.zeros(T)
    for i in range(0,len(p)):
        plot_p[i] = p[i][0]
        plot_cu[i] = cu[i][0]
    fig, axs = plt.subplots(2,1)
    axs[0].plot(range(0,T), plot_p)
    axs[0].set_title('Daily exposed')
    axs[1].plot(range(0,T), plot_cu)
    axs[1].set_title('Current cases')
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
