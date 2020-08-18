import numpy as np
from covid19_Net import covid19_Net
import covid19_R0
import graph_ops
from graph_ops import generate_random_social_graph
import datetime
from matrix_processor import matrix_processor

def repeated_runs(N, Npop, bet,tau,ph,matrices,T,sigma,filters,sim_offcampus_using=None,off_campus_mat_name="Off",sim_social=None,diagnostics=True, facility=None,quar=None, quar_filter=None, symp_test_comp=None):
    #Npop = C[0].shape[0]
    
    ns = np.zeros((N,T,9))
    cus = np.zeros((N,T,7))
    ps = np.zeros((N,T,6))
    ous = np.zeros((N,T,3))
    r0ts = np.zeros((N,T,3))
    r0infs = np.zeros((N,Npop))
    ys = np.zeros((N,Npop,1))
    xs = np.zeros((N,Npop,1))
    r0theorys = np.zeros((N,3), dtype=complex)
    norms = np.zeros((N,3), dtype=complex)
    mtss = np.zeros((N,T,1))
    
    for i in range(0,N):
        print(i)
        print(datetime.datetime.now())
        
        # To prevent updates to contact matrices from affecting future runs
        
        if sim_offcampus_using is not None:
            appt_size_rng = graph_ops.student_roommates_dist()
            Off = graph_ops.random_roommate_assignment(sim_offcampus_using,appt_size_rng)
            matrices[off_campus_mat_name] = Off
            
        if sim_social is not None:
            Soc = generate_random_social_graph(matrices['C'],sim_social['n'],sim_social['d'])
            matrices[sim_social['name']] = Soc
        
        [n,cu,p,ou,r0t,r0inf,rr,y,x] = covid19_Net(bet,tau,ph,matrices,T,sigma,filters,quar=quar,facility=facility,symp_test_comp=symp_test_comp)
        ns[i,:,:]   = n
        cus[i,:,:]  = cu
        ps[i,:,:]   = p
        ous[i,:,:]  = ou
        r0ts[i,:,:] = r0t
        r0infs[i,:] = r0inf
        ys[i,:,:]   = y
        xs[i,:,:]   = x
        
        #[n,cu,p,ou,r0t,r0inf,rr,y,x] = covid19_Net(bet,tau,ph,init,Ctmp,T,sigma,ndt,Cind)
        if diagnostics:
            proc = matrix_processor()
            Ctmp = proc.process(Npop, matrices, filters[0][0])
            [R0theory,norm,Mts,Sints,taud] = covid19_R0.covid19_R0(bet,tau,ph,Ctmp,T,sigma,r0t,n,x,rr)
            r0theorys[i,:] = R0theory
            norms[i,:] = norm
            mtss[i,:,:] = Mts
        
    return [ns,cus,ps,ous,r0ts,r0infs,ys,xs,r0theorys,norms,mtss]
    
    

def sim_with_social(N,CM_graph,RS_graph,avg_friends, dispersion, avg_contacts_per_day):
    results = []
    for x in range(N):
        results.append(generate_random_social_graph(CM_graph,RS_graph,avg_friends,dispersion,avg_contacts_per_day))
    
    return results
    
