
import matplotlib.pyplot as plt
import numpy
import scipy

# plot curves for a matrix of data
# Assumes:
#     data is a matrix.  Each row is a potential curve.  Each column is a same point in time.
#     returns a plot object.  show() will not have been called.
def plot_curves(panel, x, data, title):
    panel.set_title(title)

    n_rows = data.shape[0]
    for i in range(0,n_rows):
        panel.plot(x, data[i,:], linewidth=0.2, color="blue")
    
    
# Plot a set of data from a set of Covid simulations.
# Assumes data has been read in from scipy.io.loadmat.  This would be good to make
# more generic (Accessing state['T'] is different if gotten directly from return 
# value of simulation function than if it comes from a .mat file.  (due to some
# Matlab weirdness that doesn't allow variables to be scalars).
def gen_plots(state, prefix):
    x = range(0,state['T'][0][0])

    n = state['n']    
    fig,panels = plt.subplots(3,2, figsize=(6.5,9))
    for i in range(0,5):
        plot_curves(panels[i//2][i%2], x, n[:,:,i], "n[%d]" % i)
    fig.tight_layout()
    fig.savefig("%s_n.png" % (prefix))

    cu = state['cu']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, cu[:,:,i], "cu[%d]" % i)
    fig.tight_layout()
    fig.savefig("%s_cu.png" % (prefix))

    p = state['p']    
    fig,panels = plt.subplots(3,2, figsize=(6.5,9))
    for i in range(0,6):
        plot_curves(panels[i//2][i%2], x, p[:,:,i], "p[%d]" % i)
    fig.tight_layout()
    fig.savefig("%s_p.png" % (prefix))
    
    ou = state['ou']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, ou[:,:,i], "ou[%d]" % i)
    fig.tight_layout()
    fig.savefig("%s_ou.png" % (prefix))

    r0 = state['r0']    
    fig,panels = plt.subplots(1, figsize=(6.5,9))
    plot_curves(panels, x, r0[:,:,0], "r0[%d]" % i)
    fig.tight_layout()
    fig.savefig("%s_r0.png" % (prefix))
    
    r0t = state['r0t']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, r0t[:,:,i], "r0t[%d]" % i)
    fig.tight_layout()
    fig.savefig("%s_r0t.png" % (prefix))

    
    