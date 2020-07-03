import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import scipy.io as sio
from covid19_Net import covid19_Net
import graph_ops
from covid_simulation import repeated_runs
import datetime
import visualization as viz
import pandas as pd
import json
import sys
import stats_util as su

    
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

    class_sizes = graph_ops.calculateClassSizes(edges)

    size_limits = np.sort(np.concatenate((np.arange(10,301,10), np.arange(5,101,5))))
    #size_limits = np.array([10, 30, 50])

    trials = pd.DataFrame({
        'sz': size_limits, 
        'p_spread': 0.0,
        })
    
    su.add_summary_cols(trials,'all_final_size',0.0)
    su.add_summary_cols(trials,'spread_final_size',0.0)
    su.add_summary_cols(trials,'all_maxPrev',0.0)
    su.add_summary_cols(trials,'spread_maxPrev',0.0)
    su.add_summary_cols(trials,'all_maxSymp',0.0)
    su.add_summary_cols(trials,'spread_maxSymp',0.0)
    su.add_summary_cols(trials,'all_TimeToPeak',0.0)
    su.add_summary_cols(trials,'spread_TimeToPeak',0.0)
    su.add_summary_cols(trials,'all_TimeToNotice',0.0)
    su.add_summary_cols(trials,'spread_TimeToNotice',0.0)
    su.add_summary_cols(trials,'all_CasesAtNotice',0.0)
    su.add_summary_cols(trials,'spread_CasesAtNotice',0.0)

    print(datetime.datetime.now())

    n_students = C.shape[0]-1

    for index,row in trials.iterrows():
        include = class_sizes['Info'][class_sizes['Size'] <= row['sz']].to_list()
        edges_pruned = edges[np.logical_and(edges['Type'] == 'C',edges['Info'].isin(include))].copy()
        C2, M2, R2, S2 = graph_ops.edgesAsMatrices(edges_pruned, n_students)

        contacts = scalarC*C2 + scalarM*M + scalarR*R + scalarS*S
        
        #running the covid19_Net simulation, 30 repetitions
        C_list = [contacts,contacts] #Weekday, weekend
        C_ind = [0, 0, 0, 0, 0, 1, 1]
        [n,cu,p,ou,r0t,r0inf,y,x,r0theory,norms,mts] = repeated_runs(n_runs, bet,tau,ph,init,C_list,T,sigma,ndt,C_ind)

        max_time_step = cu.shape[1]-1
        
        fs = cu[:,max_time_step,0]
        spread = fs > 1
        n_spread = np.sum(spread)
        
        trials.at[index,'p_spread'] = np.mean(spread)

        su.calculate('all_final_size',fs)

        if n_spread > 0:
            su.calculate('spread_final_size',fs[spread])

        prevalence = p[:,:,0] + p[:,:,1] + p[:,:,2] + p[:,:,3] + p[:,:,4]
        tmp = np.max(prevalence, axis=1)

        su.calculate('all_maxPrev',tmp)

        if n_spread > 0:

            su.calculate('spread_maxPrev',tmp[spread])
        
        symp = p[:,:,1] + p[:,:,3] + p[:,:,4]
        tmp = np.max(symp, axis=1)

        su.calculate('all_maxSymp',tmp)

        if n_spread > 0:

            su.calculate('spread_maxSymp',tmp[spread])        
        
        tmp = np.argmax(prevalence, axis=1)

        su.calculate('all_TimeToPeak',tmp)

        if n_spread > 0:

            su.calculate('spread_TimeToPeak',tmp[spread])


        symp_cum = np.cumsum(symp, axis=1)
        ttn = np.sum(symp_cum==0, axis=1)
        tmp = ttn

        su.calculate('all_TimeToNotice',tmp)

        if n_spread > 0:

            su.calculate('spread_TimeToNotice',tmp[spread])
     
        basis_dates = np.clip(ttn+2, 0, max_time_step)
        cases_at_notice = prevalence[range(0,n_runs), basis_dates]
        tmp = cases_at_notice
        
        su.calculate('all_CasesAtNotice',tmp)

        if n_spread > 0:           
            su.calculate('spread_CasesAtNotice',tmp[spread])

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
        
        sio.savemat("%s_%d.mat" % (output_file,row['sz']),state)

    print(datetime.datetime.now())

    trials.to_csv("%s_summ.csv" % (output_file), index=False)
    
    # This code would load the data and generate the plots.
    # results = sio.loadmat('save_data.mat')
    # viz.gen_plots(results, "testrun")

if __name__== "__main__":
    main()

