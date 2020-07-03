
import matplotlib.pyplot as plt
import numpy
import scipy
import numpy as np

# plot curves for a matrix of data
# Assumes:
#     data is a matrix.  Each row is a potential curve.  Each column is a same point in time.
#     returns a plot object.  show() will not have been called.


def title_finder(argument,config_file): 
    switcher = { 
        "cu[0]" : config_file["cu[0]"],
        "cu[1]" : config_file["cu[1]"],
        "cu[2]" : config_file["cu[2]"],
        "n[0]" : config_file["n[0]"],
        "n[1]" : config_file["n[1]"],
        "n[2]" : config_file["n[2]"],
        "n[3]" : config_file["n[3]"],
        "n[4]" : config_file["n[4]"],
        "n[5]" : config_file["n[5]"],
        "n[6]" : config_file["n[6]"],
        "n[7]" : config_file["n[7]"],
        "n[8]" : config_file["n[8]"],
        "ou[0]" : config_file["ou[0]"],
        "ou[1]" : config_file["ou[1]"],
        "ou[2]" : config_file["ou[2]"],
        "p[0]" : config_file["p[0]"],
        "p[1]" : config_file["p[1]"],
        "p[2]" : config_file["p[2]"],
        "p[3]" : config_file["p[3]"],
        "p[4]" : config_file["p[4]"],
        "p[5]" : config_file["p[5]"],
        "r0inf" : config_file["r0inf"],
        "r0t[0]" : config_file["r0t[0]"],
        "r0t[1]" : config_file["r0t[1]"],
        "r0t[2]" : config_file["r0t[2]"], 
    } 
  
    return switcher.get(argument, argument) 


def plot_curves(panel, x, data, title,config_file):

    type =config_file['types']
    mean=config_file['mean']
    median=config_file['median']
    plot_lines_color=config_file['plot_lines_color']

    shaded_region_colors_perc90=config_file['shaded_region_colors_perc90']
    shaded_region_colors_perc50=config_file['shaded_region_colors_perc50']
    
    mean_lines_color=config_file['mean_lines_color']

    median_lines_color=config_file['median_lines_color']

    x_axis = config_file['x_axis']
    y_axis = config_file['y_axis']

    ignoreY0 = config_file['ignoreY0']


    panel.set_title(title_finder(title,config_file))

    n_rows = data.shape[0]
    
    # if title == "r0inf":
    # panel.hist(data.ravel(),range(0,51),density = True)
    # panel.set_xlabel('number of people an individual infects')
    # panel.set_ylabel('frequency(%)')

    
    
    if title == "r0inf":
        if ignoreY0:
            histData = []
            for i in range(len(data.ravel())):
                if x.ravel()[i] != 0:                   
                    histData.append(data.ravel()[i])
            panel.hist(histData,range(0,51),density = True)
            panel.set_xlabel('number of people an individual infects')
            panel.set_ylabel('frequency(%)')
            
            
        else:
            panel.hist(data.ravel(),range(0,51),density = True)
            panel.set_xlabel('number of people an individual infects')
            panel.set_ylabel('frequency(%)') 

           

    else:
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

        panel.set_xlabel(x_axis)
        panel.set_ylabel(y_axis)


def gen_plots(state, prefix,config_file):
    x = range(0,state['T'][0][0])
    
    n = state['n']    
    fig,panels = plt.subplots(3,3, figsize=(6.5,9))
    for i in range(0,9):
        plot_curves(panels[i//3][i%3], x, n[:,:,i], "n[%d]" % i,config_file)
    fig.tight_layout()
    fig.savefig("%s_n.png" % (prefix))
    

    # plot for cu
    cu = state['cu']    
    fig,panels = plt.subplots(3,2, figsize=(6.5,9))
    for i in range(0,5):
        plot_curves(panels[i//2][i%2], x, cu[:,:,i], "cu[%d]" % i,config_file)
    
    fig.tight_layout()
    fig.savefig("%s_cu.png" % (prefix))


    # cu ends

    

    p = state['p']    
    fig,panels = plt.subplots(3,2, figsize=(6.5,9))
    for i in range(0,6):
        plot_curves(panels[i//2][i%2], x, p[:,:,i], "p[%d]" % i,config_file)
    fig.tight_layout()
    fig.savefig("%s_p.png" % (prefix))
    
    ou = state['ou']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, ou[:,:,i], "ou[%d]" % i,config_file)
    fig.tight_layout()
    fig.savefig("%s_ou.png" % (prefix))


    # DC:  r0 no longer in outputs of Covid19_Net.  Revising data storage for
    # Coleman's Covid_R0 function outputs before re-adding this
    # r0 = state['r0']

    # fig,panels = plt.subplots(1, figsize=(6.5,9))
    # plot_curves(panels, x, r0[:,:,0], "r0[%d]" % i,config_file)
    # fig.tight_layout()
    # fig.savefig("%s_r0.png" % (prefix))
    


    r0t = state['r0t']    
    fig,panels = plt.subplots(3, figsize=(6.5,9))
    for i in range(0,3):
        plot_curves(panels[i], x, r0t[:,:,i], "r0t[%d]" % i,config_file)
    fig.tight_layout()
    fig.savefig("%s_r0t.png" % (prefix))

    

    r0inf = state['r0inf']
    fig,panels = plt.subplots(1, figsize=(6.5,9))  
    plot_curves(panels, state['y'], r0inf, "r0inf" ,config_file)
    fig.tight_layout()
    fig.savefig("%s_r0inf.png" % (prefix))
    
    
