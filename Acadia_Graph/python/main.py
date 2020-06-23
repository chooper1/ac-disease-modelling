import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.io as sio
from covid19_Net import covid19_Net
from graph_ops import matrixCM
from covid_simulation import repeated_runs
import datetime
import visualization as viz


    
def main():
    print(datetime.datetime.now())
    # Define parameters
    bet= [0.5,0.5,0.5]
    tau = [5, 6, 4, 8, 17, 10]
    ph = [.5, 1, 1]
    init = 1
    T=100
    sigma = [1, 1]
    ndt = 1
    n_runs = 30

    # contact matrix   
    C = matrixCM("graph_data_more.csv",12/7,1/420,500)
    print(datetime.datetime.now())
        
    #running the covid19_Net simulation, 30 repetitions
    C_list = [C,C] #Weekday, weekend
    C_ind = [0, 0, 0, 0, 0, 1, 1]
    [n,cu,p,ou,r0,r0t,y,x] = repeated_runs(n_runs, bet,tau,ph,init,C_list,T,sigma,ndt,C_ind)

    print(datetime.datetime.now())
    
    state = { 
        'bet'    : bet,
        'tau'    : tau,
        'ph'     : ph,
        'init'   : init,
        'T'      : T, 
        'sigma'  : sigma,
        'ndt'    : ndt, 
        'n_runs' : n_runs,
        'C'      : C,
        'n'      : n,
        'cu'     : cu,
        'p'      : p,
        'ou'     : ou,
        'r0'     : r0,
        'r0t'    : r0t,
        'y'      : y,
        'x'      : x
        }
        
    sio.savemat('save_data.mat',state)

    print(datetime.datetime.now())
    
    # This code would load the data and generate the plots.
    # results = sio.loadmat('save_data.mat')
    # viz.gen_plots(results, "testrun")

    return state


if __name__== "__main__":
    main()

