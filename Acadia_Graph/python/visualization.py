
import matplotlib.pyplot as plt
import numpy
import scipy
import numpy as np

# plot curves for a matrix of data
# Assumes:
#     data is a matrix.  Each row is a potential curve.  Each column is a same point in time.
#     returns a plot object.  show() will not have been called.


def plot_curves(panel, x, data, title,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color):
    panel.set_title(title)

    n_rows = data.shape[0]
    
    if type == 'lines':
        for i in range(0,n_rows):
            panel.plot(x, data[i,:], linewidth=0.2, color=plot_lines_color)

    if type == 'bands':
        
        perc90down = np.percentile(data, 5, axis=0)
        perc90up = np.percentile(data, 95, axis=0)


        panel.fill_between(x,perc90down,perc90up,color =shaded_region_colors_perc90,edgecolor = None)

        
        perc50down = np.percentile(data, 25, axis=0)
        perc50up = np.percentile(data, 75, axis=0)


        panel.fill_between(x,perc50down,perc50up,color =shaded_region_colors_perc50,edgecolor = None)
    
    if mean:
        panel.plot(x,np.mean(data,axis = 0),color = mean_lines_color)

    if median:
        panel.plot(x,np.median(data,axis = 0),color = median_lines_color)


def gen_plots(state, prefix,type ,mean ,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50 ,mean_lines_color,median_lines_color):
    x = range(0,state['T'][0][0])
    
    n = state['n']    
    fig,panels = plt.subplots(3,2, figsize=(6.5,9))
    for i in range(0,5):
        plot_curves(panels[i//2][i%2], x, n[:,:,i], "n[%d]" % i,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color)
    fig.tight_layout()
    fig.savefig("%s_n.png" % (prefix))
    

    # plot for cu
    cu = state['cu']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, cu[:,:,i], "cu[%d]" % i,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color)
    
    fig.tight_layout()
    fig.savefig("%s_cu.png" % (prefix))


    # cu ends

    

    p = state['p']    
    fig,panels = plt.subplots(3,2, figsize=(6.5,9))
    for i in range(0,6):
        plot_curves(panels[i//2][i%2], x, p[:,:,i], "p[%d]" % i,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color)
    fig.tight_layout()
    fig.savefig("%s_p.png" % (prefix))
    
    ou = state['ou']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, ou[:,:,i], "ou[%d]" % i,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color)
    fig.tight_layout()
    fig.savefig("%s_ou.png" % (prefix))

    r0 = state['r0']    
    fig,panels = plt.subplots(1, figsize=(6.5,9))
    plot_curves(panels, x, r0[:,:,0], "r0[%d]" % i,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color)
    fig.tight_layout()
    fig.savefig("%s_r0.png" % (prefix))
    
    r0t = state['r0t']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, r0t[:,:,i], "r0t[%d]" % i,type,mean,median,plot_lines_color,shaded_region_colors_perc90,shaded_region_colors_perc50,mean_lines_color,median_lines_color)
    fig.tight_layout()
    fig.savefig("%s_r0t.png" % (prefix))
    
    
