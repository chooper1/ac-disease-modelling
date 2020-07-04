import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.io as sio
from covid19_Net import covid19_Net
import graph_ops
from covid_simulation import repeated_runs
import datetime
import visualization as viz
import json
import sys


    
def main():
    edge_file   = sys.argv[1]
    output_file = sys.argv[3]

    #reading config file
    json_file = open(sys.argv[2],"r",encoding="utf-8")
    config = json.load(json_file)
    json_file.close()

    # Define parameters
    bet= config['bet']
    tau = config['tau']
    ph = config['ph']
    init = config['init']
    T=config['T']
    sigma = config['sigma']
    ndt = config['ndt']
    n_runs = config['n_runs']
    scalarC = config['scalarC']
    scalarM = config['scalarM']
    scalarR = config['scalarR']
    scalarS = config['scalarS']

    print(datetime.datetime.now())

    # load and calculate contact matrix  
    edges = graph_ops.loadEdgeList(edge_file)
    C,M,R,S = graph_ops.edgesAsMatrices(edges)
    contacts = scalarC*C + scalarM*M + scalarR*R + scalarS*S

    print(datetime.datetime.now())
        
    #running the covid19_Net simulation, 30 repetitions
    C_list = [contacts,contacts] #Weekday, weekend
    C_ind = [0, 0, 0, 0, 0, 1, 1]
    [n,cu,p,ou,r0t,r0inf,y,x,r0theory,norms,mts] = repeated_runs(n_runs, bet,tau,ph,init,C_list,T,sigma,ndt,C_ind, sim_offcampus_using=R,off_rate=2.0)

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
        'x'      : x,
        'r0theory' : r0theory,
        'norms' : norms,
        'mts' : mts,
        'mask_efficacy' : config['mask_efficacy'],
        'mask_compliance' : config['mask_compliance']
        }
        
    sio.savemat(output_file,state)

    print(datetime.datetime.now())
    
    # This code would load the data and generate the plots.
    # results = sio.loadmat('save_data.mat')
    # viz.gen_plots(results, "testrun")

if __name__== "__main__":
    main()

