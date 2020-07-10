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


def add_summary_cols(scr,name,init):
    scr[name + '_mean'] = init
    scr[name + '_min'] = init
    scr[name + '_max'] = init
    scr[name + '_median'] = init
    scr[name + '_pc5'] = init
    scr[name + '_pc25'] = init
    scr[name + '_pc75'] = init
    scr[name + '_pc95'] = init


def calculate(dataframe,index,dest,scr):
    dataframe.at[index,dest +'_mean']   = np.mean(scr)
    dataframe.at[index,dest +'_min']    = np.min(scr)
    dataframe.at[index,dest +'_max']    = np.max(scr)
    dataframe.at[index,dest +'_median'] = np.median(scr)
    dataframe.at[index,dest +'_pc5']    = np.percentile(scr, 5, axis=0)
    dataframe.at[index,dest +'_pc25']   = np.percentile(scr, 25, axis=0)
    dataframe.at[index,dest +'_pc75']   = np.percentile(scr, 75, axis=0)
    dataframe.at[index,dest +'_pc95']   = np.percentile(scr, 95, axis=0)