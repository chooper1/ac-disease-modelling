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

import sys


    
def main():
    edge_file   = sys.argv[1]
    n_runs      = int(sys.argv[2])
    output_file = sys.argv[3]

    # Define parameters
    bet= np.array([0.0080, 0.0080, 0.0080])
    tau = [5.35, 5.2, 5.0, 12.0, 12.0, 10.0]
    ph = [.5, 1, 1]
    init = 1
    T=100
    sigma = [1, 1]
    ndt = 1
    mask_efficacy = 0.44

    print(datetime.datetime.now())

    # load and calculate contact matrix  
    edges = graph_ops.loadEdgeList(edge_file)
    C,M,R,S = graph_ops.edgesAsMatrices(edges)
    contacts = 12/7*C + 1/420*M + 1/4*R + 1.0*S

    compliances = np.array([0. , 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1. ])
    
    trials = pd.DataFrame({
        'comp': compliances, 
        'p_spread': 0.0,
        'all_final_size_mean': 0.0,
        'all_final_size_min': 0.0,
        'all_final_size_max': 0.0,
        'all_final_size_median': 0.0,
        'all_final_size_pc5': 0.0,
        'all_final_size_pc25': 0.0,
        'all_final_size_pc75': 0.0,
        'all_final_size_pc95': 0.0,
        'spread_final_size_mean': 0.0,
        'spread_final_size_min': 0.0,
        'spread_final_size_max': 0.0,
        'spread_final_size_median': 0.0,
        'spread_final_size_pc5': 0.0,
        'spread_final_size_pc25': 0.0,
        'spread_final_size_pc75': 0.0,
        'spread_final_size_pc95': 0.0,
        'all_maxPrev_mean': 0.0,
        'all_maxPrev_min': 0.0,
        'all_maxPrev_max': 0.0,
        'all_maxPrev_median': 0.0,
        'all_maxPrev_pc5': 0.0,
        'all_maxPrev_pc25': 0.0,
        'all_maxPrev_pc75': 0.0,
        'all_maxPrev_pc95': 0.0,
        'spread_maxPrev_mean': 0.0,
        'spread_maxPrev_min': 0.0,
        'spread_maxPrev_max': 0.0,
        'spread_maxPrev_median': 0.0,
        'spread_maxPrev_pc5': 0.0,
        'spread_maxPrev_pc25': 0.0,
        'spread_maxPrev_pc75': 0.0,
        'spread_maxPrev_pc95': 0.0,
        'all_maxSymp_mean': 0.0,
        'all_maxSymp_min': 0.0,
        'all_maxSymp_max': 0.0,
        'all_maxSymp_median': 0.0,
        'all_maxSymp_pc5': 0.0,
        'all_maxSymp_pc25': 0.0,
        'all_maxSymp_pc75': 0.0,
        'all_maxSymp_pc95': 0.0,
        'spread_maxSymp_mean': 0.0,
        'spread_maxSymp_min': 0.0,
        'spread_maxSymp_max': 0.0,
        'spread_maxSymp_median': 0.0,
        'spread_maxSymp_pc5': 0.0,
        'spread_maxSymp_pc25': 0.0,
        'spread_maxSymp_pc75': 0.0,
        'spread_maxSymp_pc95': 0.0,
        'all_TimeToPeak_mean': 0.0,
        'all_TimeToPeak_min': 0.0,
        'all_TimeToPeak_max': 0.0,
        'all_TimeToPeak_median': 0.0,
        'all_TimeToPeak_pc5': 0.0,
        'all_TimeToPeak_pc25': 0.0,
        'all_TimeToPeak_pc75': 0.0,
        'all_TimeToPeak_pc95': 0.0,
        'spread_TimeToPeak_mean': 0.0,
        'spread_TimeToPeak_min': 0.0,
        'spread_TimeToPeak_max': 0.0,
        'spread_TimeToPeak_median': 0.0,
        'spread_TimeToPeak_pc5': 0.0,
        'spread_TimeToPeak_pc25': 0.0,
        'spread_TimeToPeak_pc75': 0.0,
        'spread_TimeToPeak_pc95': 0.0,
        'all_TimeToNotice_mean': 0.0,
        'all_TimeToNotice_min': 0.0,
        'all_TimeToNotice_max': 0.0,
        'all_TimeToNotice_median': 0.0,
        'all_TimeToNotice_pc5': 0.0,
        'all_TimeToNotice_pc25': 0.0,
        'all_TimeToNotice_pc75': 0.0,
        'all_TimeToNotice_pc95': 0.0,
        'spread_TimeToNotice_mean': 0.0,
        'spread_TimeToNotice_min': 0.0,
        'spread_TimeToNotice_max': 0.0,
        'spread_TimeToNotice_median': 0.0,
        'spread_TimeToNotice_pc5': 0.0,
        'spread_TimeToNotice_pc25': 0.0,
        'spread_TimeToNotice_pc75': 0.0,
        'spread_TimeToNotice_pc95': 0.0,
        'all_CasesAtNotice_mean': 0.0,
        'all_CasesAtNotice_min': 0.0,
        'all_CasesAtNotice_max': 0.0,
        'all_CasesAtNotice_median': 0.0,
        'all_CasesAtNotice_pc5': 0.0,
        'all_CasesAtNotice_pc25': 0.0,
        'all_CasesAtNotice_pc75': 0.0,
        'all_CasesAtNotice_pc95': 0.0,
        'spread_CasesAtNotice_mean': 0.0,
        'spread_CasesAtNotice_min': 0.0,
        'spread_CasesAtNotice_max': 0.0,
        'spread_CasesAtNotice_median': 0.0,
        'spread_CasesAtNotice_pc5': 0.0,
        'spread_CasesAtNotice_pc25': 0.0,
        'spread_CasesAtNotice_pc75': 0.0,
        'spread_CasesAtNotice_pc95': 0.0
        })
    
    print(datetime.datetime.now())
    
    appt_size_rng = graph_ops.student_roommates_dist()

    for index,row in trials.iterrows():
        bet_eff = (1 - mask_efficacy*row['comp']) * bet
        
        #running the covid19_Net simulation, 30 repetitions
        C_list = [contacts,contacts] #Weekday, weekend
        C_ind = [0, 0, 0, 0, 0, 1, 1]
        [n,cu,p,ou,r0t,r0inf,y,x,r0theory,norms,mts] = repeated_runs(n_runs, bet_eff,tau,ph,init,C_list,T,sigma,ndt,C_ind,sim_offcampus_using=R,off_rate=2.0)

        max_time_step = cu.shape[1]-1
        
        fs = cu[:,max_time_step,0]
        spread = fs > 1
        n_spread = np.sum(spread)
        
        trials.at[index,'p_spread'] = np.mean(spread)
        
        trials.at[index,'all_final_size_mean']   = np.mean(fs)
        trials.at[index,'all_final_size_min']    = np.min(fs)
        trials.at[index,'all_final_size_max']    = np.max(fs)
        trials.at[index,'all_final_size_median'] = np.median(fs)
        trials.at[index,'all_final_size_pc5']    = np.percentile(fs, 5, axis=0)
        trials.at[index,'all_final_size_pc25']   = np.percentile(fs, 25, axis=0)
        trials.at[index,'all_final_size_pc75']   = np.percentile(fs, 75, axis=0)
        trials.at[index,'all_final_size_pc95']   = np.percentile(fs, 95, axis=0)

        if n_spread > 0:
            trials.at[index,'spread_final_size_mean']   = np.mean(fs[spread])
            trials.at[index,'spread_final_size_min']    = np.min(fs[spread])
            trials.at[index,'spread_final_size_max']    = np.max(fs[spread])
            trials.at[index,'spread_final_size_median'] = np.median(fs[spread])
            trials.at[index,'spread_final_size_pc5']    = np.percentile(fs[spread], 5, axis=0)
            trials.at[index,'spread_final_size_pc25']   = np.percentile(fs[spread], 25, axis=0)
            trials.at[index,'spread_final_size_pc75']   = np.percentile(fs[spread], 75, axis=0)
            trials.at[index,'spread_final_size_pc95']   = np.percentile(fs[spread], 95, axis=0)

        prevalence = p[:,:,0] + p[:,:,1] + p[:,:,2] + p[:,:,3] + p[:,:,4]
        tmp = np.max(prevalence, axis=1)
        trials.at[index,'all_maxPrev_mean']   = np.mean(tmp)
        trials.at[index,'all_maxPrev_min']    = np.min(tmp)
        trials.at[index,'all_maxPrev_max']    = np.max(tmp)
        trials.at[index,'all_maxPrev_median'] = np.median(tmp)
        trials.at[index,'all_maxPrev_pc5']    = np.percentile(tmp, 5, axis=0)
        trials.at[index,'all_maxPrev_pc25']   = np.percentile(tmp, 25, axis=0)
        trials.at[index,'all_maxPrev_pc75']   = np.percentile(tmp, 75, axis=0)
        trials.at[index,'all_maxPrev_pc95']   = np.percentile(tmp, 95, axis=0)

        if n_spread > 0:
            trials.at[index,'spread_maxPrev_mean']   = np.mean(tmp[spread])
            trials.at[index,'spread_maxPrev_min']    = np.min(tmp[spread])
            trials.at[index,'spread_maxPrev_max']    = np.max(tmp[spread])
            trials.at[index,'spread_maxPrev_median'] = np.median(tmp[spread])
            trials.at[index,'spread_maxPrev_pc5']    = np.percentile(tmp[spread], 5, axis=0)
            trials.at[index,'spread_maxPrev_pc25']   = np.percentile(tmp[spread], 25, axis=0)
            trials.at[index,'spread_maxPrev_pc75']   = np.percentile(tmp[spread], 75, axis=0)
            trials.at[index,'spread_maxPrev_pc95']   = np.percentile(tmp[spread], 95, axis=0)
        
        symp = p[:,:,1] + p[:,:,3] + p[:,:,4]
        tmp = np.max(symp, axis=1)
        trials.at[index,'all_maxSymp_mean']   = np.mean(tmp)
        trials.at[index,'all_maxSymp_min']    = np.min(tmp)
        trials.at[index,'all_maxSymp_max']    = np.max(tmp)
        trials.at[index,'all_maxSymp_median'] = np.median(tmp)
        trials.at[index,'all_maxSymp_pc5']    = np.percentile(tmp, 5, axis=0)
        trials.at[index,'all_maxSymp_pc25']   = np.percentile(tmp, 25, axis=0)
        trials.at[index,'all_maxSymp_pc75']   = np.percentile(tmp, 75, axis=0)
        trials.at[index,'all_maxSymp_pc95']   = np.percentile(tmp, 95, axis=0)

        if n_spread > 0:
            trials.at[index,'spread_maxSymp_mean']   = np.mean(tmp[spread])
            trials.at[index,'spread_maxSymp_min']    = np.min(tmp[spread])
            trials.at[index,'spread_maxSymp_max']    = np.max(tmp[spread])
            trials.at[index,'spread_maxSymp_median'] = np.median(tmp[spread])
            trials.at[index,'spread_maxSymp_pc5']    = np.percentile(tmp[spread], 5, axis=0)
            trials.at[index,'spread_maxSymp_pc25']   = np.percentile(tmp[spread], 25, axis=0)
            trials.at[index,'spread_maxSymp_pc75']   = np.percentile(tmp[spread], 75, axis=0)
            trials.at[index,'spread_maxSymp_pc95']   = np.percentile(tmp[spread], 95, axis=0)        
        
        tmp = np.argmax(prevalence, axis=1)
        trials.at[index,'all_TimeToPeak_mean']   = np.mean(tmp)
        trials.at[index,'all_TimeToPeak_min']    = np.min(tmp)
        trials.at[index,'all_TimeToPeak_max']    = np.max(tmp)
        trials.at[index,'all_TimeToPeak_median'] = np.median(tmp)
        trials.at[index,'all_TimeToPeak_pc5']    = np.percentile(tmp, 5, axis=0)
        trials.at[index,'all_TimeToPeak_pc25']   = np.percentile(tmp, 25, axis=0)
        trials.at[index,'all_TimeToPeak_pc75']   = np.percentile(tmp, 75, axis=0)
        trials.at[index,'all_TimeToPeak_pc95']   = np.percentile(tmp, 95, axis=0)

        if n_spread > 0:
            trials.at[index,'spread_TimeToPeak_mean']   = np.mean(tmp[spread])
            trials.at[index,'spread_TimeToPeak_min']    = np.min(tmp[spread])
            trials.at[index,'spread_TimeToPeak_max']    = np.max(tmp[spread])
            trials.at[index,'spread_TimeToPeak_median'] = np.median(tmp[spread])
            trials.at[index,'spread_TimeToPeak_pc5']    = np.percentile(tmp[spread], 5, axis=0)
            trials.at[index,'spread_TimeToPeak_pc25']   = np.percentile(tmp[spread], 25, axis=0)
            trials.at[index,'spread_TimeToPeak_pc75']   = np.percentile(tmp[spread], 75, axis=0)
            trials.at[index,'spread_TimeToPeak_pc95']   = np.percentile(tmp[spread], 95, axis=0)


        symp_cum = np.cumsum(symp, axis=1)
        ttn = np.sum(symp_cum==0, axis=1)
        tmp = ttn
        trials.at[index,'all_TimeToNotice_mean']   = np.mean(tmp)
        trials.at[index,'all_TimeToNotice_min']    = np.min(tmp)
        trials.at[index,'all_TimeToNotice_max']    = np.max(tmp)
        trials.at[index,'all_TimeToNotice_median'] = np.median(tmp)
        trials.at[index,'all_TimeToNotice_pc5']    = np.percentile(tmp, 5, axis=0)
        trials.at[index,'all_TimeToNotice_pc25']   = np.percentile(tmp, 25, axis=0)
        trials.at[index,'all_TimeToNotice_pc75']   = np.percentile(tmp, 75, axis=0)
        trials.at[index,'all_TimeToNotice_pc95']   = np.percentile(tmp, 95, axis=0)

        if n_spread > 0:
            trials.at[index,'spread_TimeToNotice_mean']   = np.mean(tmp[spread])
            trials.at[index,'spread_TimeToNotice_min']    = np.min(tmp[spread])
            trials.at[index,'spread_TimeToNotice_max']    = np.max(tmp[spread])
            trials.at[index,'spread_TimeToNotice_median'] = np.median(tmp[spread])
            trials.at[index,'spread_TimeToNotice_pc5']    = np.percentile(tmp[spread], 5, axis=0)
            trials.at[index,'spread_TimeToNotice_pc25']   = np.percentile(tmp[spread], 25, axis=0)
            trials.at[index,'spread_TimeToNotice_pc75']   = np.percentile(tmp[spread], 75, axis=0)
            trials.at[index,'spread_TimeToNotice_pc95']   = np.percentile(tmp[spread], 95, axis=0)

      
        basis_dates = np.clip(ttn+2, 0, max_time_step)
        cases_at_notice = prevalence[range(0,n_runs), basis_dates]
        tmp = cases_at_notice
        trials.at[index,'all_CasesAtNotice_mean']   = np.mean(tmp)
        trials.at[index,'all_CasesAtNotice_min']    = np.min(tmp)
        trials.at[index,'all_CasesAtNotice_max']    = np.max(tmp)
        trials.at[index,'all_CasesAtNotice_median'] = np.median(tmp)
        trials.at[index,'all_CasesAtNotice_pc5']    = np.percentile(tmp, 5, axis=0)
        trials.at[index,'all_CasesAtNotice_pc25']   = np.percentile(tmp, 25, axis=0)
        trials.at[index,'all_CasesAtNotice_pc75']   = np.percentile(tmp, 75, axis=0)
        trials.at[index,'all_CasesAtNotice_pc95']   = np.percentile(tmp, 95, axis=0)

        if n_spread > 0:
            trials.at[index,'spread_CasesAtNotice_mean']   = np.mean(tmp[spread])
            trials.at[index,'spread_CasesAtNotice_min']    = np.min(tmp[spread])
            trials.at[index,'spread_CasesAtNotice_max']    = np.max(tmp[spread])
            trials.at[index,'spread_CasesAtNotice_median'] = np.median(tmp[spread])
            trials.at[index,'spread_CasesAtNotice_pc5']    = np.percentile(tmp[spread], 5, axis=0)
            trials.at[index,'spread_CasesAtNotice_pc25']   = np.percentile(tmp[spread], 25, axis=0)
            trials.at[index,'spread_CasesAtNotice_pc75']   = np.percentile(tmp[spread], 75, axis=0)
            trials.at[index,'spread_CasesAtNotice_pc95']   = np.percentile(tmp[spread], 95, axis=0)
        

                


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
            'mts' : mts
            }
        
        sio.savemat("%s_%d.mat" % (output_file,row['sz']),state)

    print(datetime.datetime.now())

    trials.to_csv("%s_summ.csv" % (output_file), index=False)
    
    # This code would load the data and generate the plots.
    # results = sio.loadmat('save_data.mat')
    # viz.gen_plots(results, "testrun")

if __name__== "__main__":
    main()

