


import pytools as pt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path
import copy

import os, sys
from read_vsc_data import *    # read_vsc_data reads txt file output from vlsvintpol.py into a dictionary




# make a function that takes a time series (time, data) and the time of the pulse transition and calculates before and after values. There may be a couple different ways to do this

def timeseries_before_after(t, data, t_mid):
    ''' t is array of times
        data is the timeseries, same size as t
        t_mid is the time of the pulse transition
        Returns a 2-element tuple: (before, after) the nominal values of input data at the beginning and end of the pulse
    '''
    mask_before = (t <= t_mid)
    mask_after = (t > t_mid)
    mean_before = np.nanmean(data[mask_before])
    mean_after = np.nanmean(data[mask_after])
    std_before = np.nanstd(data[mask_before])
    std_after = np.nanstd(data[mask_after])
    return mean_before, mean_after, std_before, std_after


# MAIN:

# Vlasiator global run data

f_EGL = pt.vlsvfile.VlsvReader( get_vlsvfile_fullpath("EGL", 1000) )
f_EGP = pt.vlsvfile.VlsvReader( get_vlsvfile_fullpath("EGP", 300) )

run_velocity = {"EGL":-750000,"EGP":-1000000}
run_t0 = {"EGL":620, "EGP":269 }
run_xmax = {"EGL":f_EGL.read_parameter('xmax'), "EGP":f_EGP.read_parameter('xmax')}


# LOAD Vlasiator data

run = 'EGL'    # 'EGP', ...
root_dir = '/wrk-vakka/users/horakons/carrington/plots/'
save_dir = '{}{}/timeseries/'.format(root_dir, run.upper())         #CUSTOM

dct_vsc, nvars_vsc = read_vsc_data('txt_files/virtual_sc_data_vlasiator_comp_{}.txt'.format(run.lower()))

velocity = run_velocity[run]
t0 = run_t0[run]
xmax = run_xmax[run]

R_EARTH = 6.371e6            #check what value is used in simulations

# every item is a 2d array [x,y] where x indexes time, y indexes position

t_vsc = dct_vsc['t'][:,0]
X_RE = dct_vsc['X_RE'][0,:]
Y_RE = dct_vsc['Y_RE'][0,:]
Z_RE = dct_vsc['Z_RE'][0,:]
CELLID = dct_vsc['CELLID'][0,:]
del dct_vsc['t']
del dct_vsc['X_RE']
del dct_vsc['Y_RE']
del dct_vsc['Z_RE']
del dct_vsc['CELLID']


ntime = t_vsc.size
npos_vsc = X_RE.size
#ntime = t_vsc.shape[0]
#npos_vsc = t_vsc.shape[1]

# LOAD pulse data (sc)

df_sc = pd.read_csv('csv_files/pulse_omni.csv')
df_sc = df_sc.sort_values(by=['PULSE', 'TIME_PULSE'])
npts_per_pulse = np.array(np.where(df_sc['PULSE'] == 0)).size
n_pulse = int(df_sc.shape[0] / npts_per_pulse)
dct_sc = {}
for varname in df_sc:
    arr = np.array(df_sc[varname])
    dct_sc[varname] = [arr[i:i + npts_per_pulse] for i in range(0, len(arr), npts_per_pulse)]   #list comprehension: makes a list of pulses with equal length
    #exec( 'arr = np.array(df_sc[varname])' )
    #exec( '{}=[arr[i:i + npts_per_pulse] for i in range(0, len(arr), npts_per_pulse)]'.format(varname) )
    # This creates a list for each varname, each element of the list represents a different pulse and has pulse_len elements
    # example: dct_sc['PULSE'][0] = [0, 0, 0, 0, 0,... 0], dct_sc['PULSE'][1] = [1, 1, 1, 1, ... 1]
    # dct_sc['TIME_PULSE'] is a list of times for the pulses (in seconds), etc.


pulse_len_sec = dct_sc['TIME_PULSE'][0][-1] - dct_sc['TIME_PULSE'][0][0]

# dct keys:
# 'B.x', 'B.y', 'B.z', 'proton/rho', 'proton/temperature', 'proton/vg_beta', 'proton/vg_pdyn', 'proton/vg_v.magnitude', 'proton/vg_v.x', 'proton/vg_v.y', 'proton/vg_v.z', 'vg_b_vol.magnitude', 'vg_e_vol.magnitude'
# dct_sc names:
# PDYN  FLOW_SPEED        T  BZ_GSE        B   TIME_PULSE  PULSE

plotvar_sc = ['PDYN', 'FLOW_SPEED', 'T', 'BZ_GSE', 'B']
plotvar_vsc = ['proton/vg_pdyn', 'proton/vg_v.magnitude', 'proton/temperature', 'B.z', 'vg_b_vol.magnitude']

# TODO:
# *save t_mid and t0 into lists ... change their names too
# *make a function that takes a time series (time, data) and the time of the pulse transition and calculates before and after values. There may be a couple different ways to do this
# *using the above function calculate before and afters for all variables, for pulses in the SC data and for all positions in VSC data
# Designate a virtual spacecraft position to be used as a reference
# write a function that takes the matching variable data and computes the 'best' match (normalize to 1?), either based on before/after or based on Pearson coefficient or something like that
# write a function that computes the 'best' match overall across multiple variables. Specify a virtual spacecraft node

# plot only the time series that best matches overall (or the top 2-3 best matches)

dct_before_after_vsc = { 't_mid_list': [], 
                         'mean_before_list': [],  # (npos_vsc)-element list of arrays with len(plotvar_vsc) elements
                         'mean_after_list': [],   # ''
                         'std_before_list': [],   # ''
                         'std_after_list': [],    # ''
                         }

dct_before_after_sc = copy.deepcopy(dct_before_after_vsc)

#t_mid_vsc_list = []
#t_mid_sc_list = []         # convert from epoch (msec) to SI (sec)
#data_before_sc_list = []  # (n_pulse)-element list of arrays each with len(plotvar_sc) elements
#data_after_sc_list = []  # (n_pulse)-element list of arrays each with len(plotvar_sc) elements
#data_before_vsc_list = []  # (npos_vsc)-element list of arrays with len(plotvar_vsc) elements (which is same as _sc_list above)
#data_after_vsc_list = []  # (npos_vsc)-element list of arrays with len(plotvar_vsc) elements (which is same as _sc_list above)


for i in range(npos_vsc):
    t_mid_temp = t0 + (X_RE[i]*R_EARTH - xmax) / velocity
    #t_mid_temp = t0 + (X_RE[0,i]*R_EARTH - xmax) / velocity
    dct_before_after_vsc['t_mid_list'].append( t_mid_temp )
    dct_before_after_vsc['mean_before_list'].append(np.ndarray(len(plotvar_vsc)))
    dct_before_after_vsc['mean_after_list'].append(np.ndarray(len(plotvar_vsc)))
    dct_before_after_vsc['std_before_list'].append(np.ndarray(len(plotvar_vsc)))
    dct_before_after_vsc['std_after_list'].append(np.ndarray(len(plotvar_vsc)))
    for j in range(len(plotvar_vsc)):
        mean_before, mean_after, std_before, std_after = timeseries_before_after(t_vsc, dct_vsc[plotvar_vsc[j]][:,i], t_mid_temp)
        dct_before_after_vsc['mean_before_list'][i][j] = mean_before
        dct_before_after_vsc['mean_after_list'][i][j] = mean_after
        dct_before_after_vsc['std_before_list'][i][j] = std_before
        dct_before_after_vsc['std_after_list'][i][j] = std_after


for i in range(n_pulse):
    t_mid_temp = np.average(dct_sc['TIME_PULSE'][i])
    dct_before_after_sc['t_mid_list'].append( t_mid_temp )
    dct_before_after_sc['mean_before_list'].append(np.ndarray(len(plotvar_sc)))
    dct_before_after_sc['mean_after_list'].append(np.ndarray(len(plotvar_sc)))
    dct_before_after_sc['std_before_list'].append(np.ndarray(len(plotvar_sc)))
    dct_before_after_sc['std_after_list'].append(np.ndarray(len(plotvar_sc)))
    for j in range(len(plotvar_sc)):
        mean_before, mean_after, std_before, std_after = timeseries_before_after(dct_sc['TIME_PULSE'][i], dct_sc[plotvar_sc[j]][i], t_mid_temp)
        dct_before_after_sc['mean_before_list'][i][j] = mean_before
        dct_before_after_sc['mean_after_list'][i][j] = mean_after
        dct_before_after_sc['std_before_list'][i][j] = std_before
        dct_before_after_sc['std_after_list'][i][j] = std_after



# Designate a virtual spacecraft position to be used as a reference

#vsc_index = 5            # reimpliment this later... find the nearest point to specified coordinates
vsc_index = 0            # reimpliment this later... find the nearest point to specified coordinates


chisq_dof = np.ndarray(n_pulse)

dof = len(plotvar_sc) * 2    # *2 comes from before and after parts
for i in range(n_pulse):
    diff_after = dct_before_after_sc['mean_after_list'][i] - dct_before_after_vsc['mean_after_list'][vsc_index]
    diff_before = dct_before_after_sc['mean_before_list'][i] - dct_before_after_vsc['mean_before_list'][vsc_index]
    std_after = (dct_before_after_sc['std_after_list'][i]**2 + dct_before_after_vsc['std_after_list'][vsc_index]**2 )**0.5
    std_before = (dct_before_after_sc['std_before_list'][i]**2 + dct_before_after_vsc['std_before_list'][vsc_index]**2 )**0.5
    chisq_dof[i] = np.sum((diff_after / std_after)**2 + (diff_before / std_before)**2 ) / dof



ind_sc_best_all = np.argsort(chisq_dof)
nbest = 1   # number of plots to make
ind_sc_best = ind_sc_best_all[0:nbest]


#ind_sc_bestmatch = np.where(chisq_dof == np.nanmin(chisq_dof))[0][0]






#PLOT timeseries for and Vlasiator and the best-matching SC pulse
# make a time-series plot at the given coordinates of all the variables
# test for vector components like B.x, B.y, B.z in a simple way: keep track of the last plotted thing,
# if it was NAME.x and this is NAME.y (look for the '.') then don't close the current plot.
# add a legend regardless of whether multiple things were plotted
#keys = sorted(dct_vsc.keys())    # sort them alphabetically
for i in range(npos_vsc):
    t_mid_vsc = dct_before_after_vsc['t_mid_list'][i]
    mask_vsc = ((t_vsc >= t_mid_vsc-pulse_len_sec/2) & (t_vsc < t_mid_vsc + pulse_len_sec/2))
    nvar = len(plotvar_sc)
    fig, axes = plt.subplots(nvar, 1)      #axes is a 1D array with npos_vsc elements
    for j in range(nvar):
        key_sc = plotvar_sc[j]
        key_vsc = plotvar_vsc[j]
        value_vsc = dct_vsc[key_vsc]
        for ind_sc in ind_sc_best:
            axes[j].plot( (dct_sc['TIME_PULSE'][ind_sc] - dct_before_after_sc['t_mid_list'][ind_sc] ), dct_sc[key_sc][ind_sc], color = 'green')
        axes[j].plot(t_vsc[mask_vsc]-t_mid_vsc, value_vsc[mask_vsc, i], color = 'red', label = key_vsc.replace('_', '\\_'))
        minv = np.nanmin(dct_sc[key_sc][ind_sc])
        axes[j].set_ylim([np.nanmin([0, minv]), None])
        axes[j].legend()
        if j < nvar-1:
            for xlabel_j in axes[j].get_xticklabels():
                xlabel_j.set_visible(False) 
    axes[j].xlabel = 'time [seconds]'
    plt.subplots_adjust(hspace=0)
    filename = '{}sc_vlasiator_comp_bestmatch_{}_x{}y{}z{}.png'.format(save_dir, run, X_RE[i], Y_RE[i], Z_RE[i])
    mkdir_path(filename)
    plt.savefig(filename)
    plt.close()




#PLOT timeseries for ALL pulses, and overplot the Vlasiator results for the same variables
for i in range(npos_vsc):
    # make a time-series plot at the given coordinates of all the variables
    # test for vector components like B.x, B.y, B.z in a simple way: keep track of the last plotted thing,
    # if it was NAME.x and this is NAME.y (look for the '.') then don't close the current plot.
    # add a legend regardless of whether multiple things were plotted
    #keys = sorted(dct_vsc.keys())    # sort them alphabetically
    t_mid_vsc = dct_before_after_vsc['t_mid_list'][i]
    #mask_vsc = ((t_vsc[:,0] >= t_mid_vsc-pulse_len_sec/2) & (t_vsc[:,0] < t_mid_vsc + pulse_len_sec/2))
    mask_vsc = ((t_vsc >= t_mid_vsc-pulse_len_sec/2) & (t_vsc < t_mid_vsc + pulse_len_sec/2))
    nvar = len(plotvar_sc)
    fig, axes = plt.subplots(nvar, 1)      #axes is a 1D array with npos_vsc elements
    for j in range(nvar):
        key_sc = plotvar_sc[j]
        key_vsc = plotvar_vsc[j]
        value_vsc = dct_vsc[key_vsc]
        for k in range(n_pulse):
            axes[j].plot( (dct_sc['TIME_PULSE'][k] - dct_before_after_sc['t_mid_list'][k] ), dct_sc[key_sc][k], color = 'green')
            #exec( "axes[j].plot( (TIME_PULSE[k] - t_mid_sc_list[k] ), {}[k], color = 'green')".format(key_sc) )
        #axes[j].plot(t[:,i], value_vsc[:,i])
        axes[j].plot(t_vsc[mask_vsc]-t_mid_vsc, value_vsc[mask_vsc,i], color = 'red', label = key_vsc.replace('_', '\\_'))
        #axes[j].plot(t_vsc[mask_vsc,i]-t_mid_vsc, value_vsc[mask_vsc,i], color = 'red', label = key_vsc.replace('_', '\\_'))
        axes[j].set_ylim([0, None])
        axes[j].legend()
        if j < nvar-1:
            for xlabel_j in axes[j].get_xticklabels():
                xlabel_j.set_visible(False)
    axes[j].xlabel = 'time [seconds]'
    plt.subplots_adjust(hspace=0)
    #fig.tight_layout()
    filename = '{}sc_vlasiator_comp_{}_x{}y{}z{}.png'.format(save_dir, run, X_RE[i], Y_RE[i], Z_RE[i])
    #filename = '{}sc_vlasiator_comp_{}_x{}y{}z{}.png'.format(save_dir, run, X_RE[0,i], Y_RE[0,i], Z_RE[0,i])
    mkdir_path(filename)
    plt.savefig(filename)
    plt.close()




