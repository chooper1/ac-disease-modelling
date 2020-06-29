import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.io as sio
from covid19_Net import covid19_Net
import graph_ops
from covid_simulation import repeated_runs
import datetime
import visualization as viz

import sys


    
def main():
    edge_file   = sys.argv[1]
    n_runs      = int(sys.argv[2])
    output_file = sys.argv[3]

    # Define parameters
    bet= [0.5,0.5,0.5]
    tau = [5, 6, 4, 8, 17, 10]
    ph = [.5, 1, 1]
    init = 1
    T=100
    sigma = [1, 1]
    ndt = 1

    print(datetime.datetime.now())

    # load and calculate contact matrix  
    edges = graph_ops.loadEdgeList(edge_file)
    C,M,R,S = graph_ops.edgesAsMatrices(edges)
    contacts = 12/7*C + 1/420*M + 1/4*R + 1*S

    print(datetime.datetime.now())
        
    #running the covid19_Net simulation, 30 repetitions
    C_list = [contacts,contacts] #Weekday, weekend
    C_ind = [0, 0, 0, 0, 0, 1, 1]
    [n,cu,p,ou,r0t,r0inf,y,x] = repeated_runs(n_runs, bet,tau,ph,init,C_list,T,sigma,ndt,C_ind)

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
        'C'      : contacts,
        'n'      : n,
        'cu'     : cu,
        'p'      : p,
        'ou'     : ou,
        'r0t'    : r0t,
        'r0inf'  : r0inf,
        'y'      : y,
        'x'      : x
        }
        
    sio.savemat(output_file,state)

    print(datetime.datetime.now())
    
    # This code would load the data and generate the plots.
    # results = sio.loadmat('save_data.mat')
    # viz.gen_plots(results, "testrun")

if __name__== "__main__":
    main()

