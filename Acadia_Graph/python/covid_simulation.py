import numpy as np
from covid19_Net import covid19_Net
from graph_ops import generate_random_social_graph
import datetime

def repeated_runs(N, bet,tau,ph,init,C,T,sigma,ndt,nmatrices):
    Npop = C[0].shape[0]
    
    ns = np.zeros((N,T,5))
    cus = np.zeros((N,T,3))
    ps = np.zeros((N,T,6))
    ous = np.zeros((N,T,3))
    r0s = np.zeros((N,T,1))
    r0ts = np.zeros((N,T,3))
    ys = np.zeros((N,Npop,1))
    xs = np.zeros((N,Npop,1))
    
    for i in range(0,N):
        print(i)
        print(datetime.datetime.now())
        [n,cu,p,ou,r0,r0t,y,x] = covid19_Net(bet,tau,ph,init,C,T,sigma,ndt,nmatrices)
        print(datetime.datetime.now())
        ns[i,:,:]   = n
        cus[i,:,:]  = cu
        ps[i,:,:]   = p
        ous[i,:,:]  = ou
        r0s[i,:,:]  = r0
        r0ts[i,:,:] = r0t
        ys[i,:,:]   = y
        xs[i,:,:]   = x
        
    return [ns,cus,ps,ous,r0s,r0ts,ys,xs]
    
    

def sim_with_social(N,CM_graph,RS_graph,avg_friends, dispersion, avg_contacts_per_day):
    results = []
    for x in range(N):
        results.append(generate_random_social_graph(CM_graph,RS_graph,avg_friends,dispersion,avg_contacts_per_day))
    
    return results
    
