#this file provides an example for running the code to get the number of people infected
#before the first symptomatic case

import matplotlib.pyplot as plt
import numpy as np
from covid19_Net_First_Symptomatic import covid19_Net_First_Symptomatic

def main():
    # Define parameters
    bet= [0.5,0.5,0.5]
    tau = [5, 6, 4, 8, 17, 10]
    ph = [.5, 1, 1]
    init = 1
    T=100
    Npop = 4000
    sigma = [1, 1]
    ndt = 1

    # contact matrix
    C_Weekday=10*np.random.uniform(0,1,size=(Npop,Npop)) / Npop
    C_Weekend=10*np.random.uniform(0,1,size=(Npop,Npop)) / Npop

    C = [C_Weekday, C_Weekend]
    C_ind = [0, 0, 0, 0, 0, 1, 1]
    #running the covid19_Net simulation
    [n,cu,p,ou,r0,r0t,y,x,LastTime] = covid19_Net_First_Symptomatic(bet,tau,ph,init,C,T,sigma,ndt,C_ind)


    #plotting
    plot_n = np.zeros(LastTime)
    plot_cu = np.zeros(LastTime)
    for i in range(0,LastTime):
        plot_n[i] = n[i][0]
        plot_cu[i] = cu[i][0]

    fig, axs = plt.subplots(2,1)
    axs[0].plot(range(0,LastTime), plot_n)
    axs[0].set_title('New Cases')
    axs[1].plot(range(0,LastTime), plot_cu)
    axs[1].set_title('Current cases')
    plt.show()

    #cu[:][0] gives the total number of cases before the
    #first symptomatic case appeared

if __name__== "__main__":
    main()
