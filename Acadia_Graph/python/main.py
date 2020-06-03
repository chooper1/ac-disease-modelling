import matplotlib.pyplot as plt
import numpy as np
from covid19_Net import covid19_Net

def main():
    # Define parameters
    bet= [0.5,0.5,0.5]
    tau = [5, 6, 4, 8, 17, 10]
    ph = [.75, 1, 0.9]
    init = 4
    T=100
    Npop = 1000
    sigma = [1, 1]
    ndt = 1
    # contact matrix
    C=10*np.random.uniform(0,1,size=(Npop,Npop)) / Npop

    #running the covid19_Net simulation
    [n,cu,p,ou,r0,r0t,y,x] = covid19_Net(bet,tau,ph,init,C,T,sigma,ndt)

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


if __name__== "__main__":
    main()
