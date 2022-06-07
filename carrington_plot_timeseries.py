import pytools as pt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy.interpolate import griddata
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path

import os, sys
from read_vsc_data import *    # read_vsc_data reads txt file output from vlsvintpol.py into a dictionary
import pytools as pt
import numpy as np

import sys
import argparse
import keogram as keogram
from copy import deepcopy



# Program purpose: makes timeseries plots of the main inferred parameters indicative of solar storms

parser = argparse.ArgumentParser()
parser.add_argument('-run', default='EGL', help="the Vlasiator run, e.g. egl or EGL" )
parser.add_argument('-n_sc', help="number of coordinates in the data files" )
parser.add_argument('-t', default='t_sec', help="number of coordinates in the data files" )
args = parser.parse_args()

#import pdb; pdb.set_trace()

root_dir = '/wrk-vakka/users/horakons/carrington/data/'
save_dir = '{}{}/aurora_plot_data/'.format(root_dir, args.run.upper())         #CUSTOM

datafile = save_dir + 'data.csv'

plot_dir = '/wrk-vakka/users/horakons/carrington/plots/' + args.run.upper() + '/summary/'
plot_file = plot_dir  + 'carrington_summary_timeseries_{}.png'.format(args.run) 


if args.n_sc is None:
    n_sc_input = None
else:
    n_sc_input= int(args.n_sc)


dct, nvars = read_vsc_data(datafile, n_sc = n_sc_input, tvar = args.t)


t = dct[args.t][:,0]
ntimes = t.size


# TODO: implement this separately for North and South poles

proton_DEF_15811_eV_max = np.zeros([ntimes])
last_closed_north_deg = np.zeros([ntimes])
last_closed_south_deg = np.zeros([ntimes])
first_open_north_deg = np.zeros([ntimes])
first_open_south_deg = np.zeros([ntimes])
jpar_north_max = np.zeros([ntimes])
jpar_south_max = np.zeros([ntimes])
magnetopause_standoff_RE = np.zeros([ntimes])

# calculate quantities of interest

theta_GSE_deg = dct['theta_GSE_deg'][0,:]
phi_GSE_deg = dct['phi_GSE_deg'][0,:]

mask_north = theta_GSE_deg < 90
mask_south = theta_GSE_deg > 90


#theta_GSE_deg goes from 0 (north pole)--> 180 (south pole)
#phi_GSE_deg goes from -180 to 180  (0 = noon)


# Key names:
#r_GSE_m,theta_GSE_deg,phi_GSE_deg,r_GSE_traced_m,theta_GSE_traced_deg,phi_GSE_traced_deg,t_sec,open_vs_closed,proton_DEF_00500_eV,proton_DEF_00889_eV,proton_DEF_01581_eV,proton_DEF_02811_eV,proton_DEF_05000_eV,proton_DEF_08891_eV,proton_DEF_15811_eV,proton_DEF_28117_eV,proton_DEF_50000_eV,$(B_0 / B)J_\parallel$  $[A/km^2]$,Magnetopause_radius_RE

keys = ['proton_DEF_15811_eV', 'open_vs_closed', '$(B_0 / B)J_\parallel$  $[A/km^2]$', 'Magnetopause_radius_RE\n']


for i in range(ntimes):
    # 1. proton DEF
    proton_DEF_15811_eV_max[i] = np.nanmax(dct[keys[0]][i,:])
    # 2. last closed field line
    #inds_closed_north, = np.where( (dct[keys[1]][i,:] == 0.0) & (theta_GSE_deg < 90))
    #inds_closed_south, = np.where( (dct[keys[1]][i,:] == 0.0) & (theta_GSE_deg > 90))
    #last_closed_north_deg[i] = 90-np.nanmin(theta_GSE_deg[inds_closed_north])
    #last_closed_south_deg[i] = 90-np.nanmax(theta_GSE_deg[inds_closed_south]) 
    inds_open_north, = np.where( (dct[keys[1]][i,:] == 1.0) & (theta_GSE_deg < 90))
    inds_open_south, = np.where( (dct[keys[1]][i,:] == 1.0) & (theta_GSE_deg > 90))
    first_open_north_deg[i] = 90-np.nanmax(theta_GSE_deg[inds_open_north])
    first_open_south_deg[i] = 90-np.nanmin(theta_GSE_deg[inds_open_south]) 
    # 3. parallel current
    jpar_north_max[i] = np.nanmax( np.abs(dct[keys[2]][i,mask_north]) )
    jpar_south_max[i] = np.nanmax( np.abs(dct[keys[2]][i,mask_south]) )
    # 4. magnetopause standoff
    magnetopause_standoff_RE[i] = np.nanmin( np.abs(dct[keys[3]][i,:]) )


dct_plot = {}
dct_plot['Proton DEF, 15811 eV (max.)'] = proton_DEF_15811_eV_max
#dct_plot['Last closed lat. (deg.)'] = last_closed_north_deg
dct_plot['First open lat. (deg.)'] = first_open_north_deg
dct_plot['|J| (max.)'] = jpar_north_max
dct_plot['Magnetopause standoff (RE)'] = magnetopause_standoff_RE



dct_plot_keys = dct_plot.keys()
nvar = len(dct_plot_keys)

fig, axes = plt.subplots(nvar, 1)      #axes is a 1D array with npos elements

for i, key in enumerate(dct_plot_keys):
    label = key.replace('_', '\\_')
    axes[i].plot(t, dct_plot[key], label = label)
    axes[i].legend()
    if i < nvar-1:
        for xlabel_i in axes[i].get_xticklabels():
            xlabel_i.set_visible(False)
    axes[i].xlabel = 'time [seconds]'
    plt.subplots_adjust(hspace=0)
    #fig.tight_layout()
print(plot_file)
mkdir_path(plot_file)
plt.savefig(plot_file)
plt.close()

save(save_dir+'carrington_plot_timeseries.pickle', t=t, dct=dct_plot)












