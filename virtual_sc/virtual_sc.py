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
from time import time

from scipy.signal import savgol_filter


def squareAxes(r1, r2):
    # assume r1 and r2 are 1D arrays of the same length...
    # ...but they contain the coordinates of a 2D array
    # return the same arrays, reshaped as 2D arrays
    arr1 = np.unique(r1)
    arr2 = np.unique(r2)
    n1=arr1.size
    n2=arr2.size
    if n1*n2 != len(r1):
        print("Error: coordinates r1 and r2 aren't arranged in a regular grid!")
        return None, None
    else:
        if r1[1] == r1[0]:
            dim_x = n1
            dim_y = n2
        else:
            dim_x = n2
            dim_y = n1
        r1_2d = r1.reshape(dim_x, dim_y)
        r2_2d = r2.reshape(dim_x, dim_y)
        return r1_2d, r2_2d



def make_label(key, filter=False, dt=False, abs=False, norm=False):
    label = key.replace('_', '\\_')
    if filter:
        label = r"$\Delta$ " + label
    if dt:
        label = r"$\partial_t$ " + label
    if abs:
        label = 'modulus ' + label
    if norm:
        label = label + ' (norm)'
    return label


def str2float(string):
    # convert string to float
    # need to handle the case where input string is None
    if string is None:
        return None
    else:
        return float(string)


#Parse command line arguments:

parser = argparse.ArgumentParser()
parser.add_argument('-run', default='EGL', help="the Vlasiator run, e.g. egl or EGL" )
parser.add_argument('-nproc', default=1, help="number of processes, multithreading" )
parser.add_argument('-filein', help="name of the input file" )
parser.add_argument('-fileprefix', default='timeseries', help="prefix of the output png file" )
parser.add_argument('-n_sc', help="(if necessary) specify the number of virtual spacecraft " )
parser.add_argument('-t', default='t', help="time variable name" )
parser.add_argument('-r1', default='X_RE', help="first spatial coordinate, variable name" )
parser.add_argument('-r2', default='Y_RE', help="second spatial coordinate, variable name" )
parser.add_argument('-r3', default='Z_RE', help="third spatial coordinate, variable name" )
parser.add_argument('-var', nargs='*', help="list of variable names to plot")

#parser.add_argument('-xvar', help="x-axis variable of the keogram, as listed in file (default: integrate position vector)" )
#parser.add_argument('-yvar', default='t', help="y-axis variable of the keogram, as listed in file (default: t)" )
parser.add_argument('-dl', action='store_true', help="if set, compute curvilinear distances (for keogram), assuming r1,r2,r3 are cartesian coordinates ")
parser.add_argument('-ta', action='store_true', help="if set, compute alfven crossing time. Assume coordinates form a curvilinear path.")
parser.add_argument('-dt', action='store_true', help="if set, compute the time derivative of the variable")
parser.add_argument('-abs', action='store_true', help="if set, compute the absolute value of the final result")
parser.add_argument('-xlim', nargs='*', help="axes xlim (list [xmin, xmax])" )
parser.add_argument('-ylim', nargs='*', help="axes ylim (list [xmin, xmax])" )
#parser.add_argument('-vmin', nargs='*', help="list of mins (one for each variable)" )
#parser.add_argument('-vmax', nargs='*', help="list of maxes (one for each variable)" )
parser.add_argument('-vmin', help="min value of variable to plot (use case: when -var specifies only one variable)" )
parser.add_argument('-vmax', help="max value of variable to plot (use case: when -var specifies only one variable)" )
parser.add_argument('-filter', action='store_true',  help="apply Savitsky-Golay filter to all data, subtract filtered data to get delta")
parser.add_argument('-log', action='store_true',  help="plot logarithm of z-data")
parser.add_argument('-norm', action='store_true',  help="if filter is set, divide filtered data by trend. if filter is not set, divide by background (initial condition)")
parser.add_argument('-timeseries', action='store_true',  help="make time series plots of all variables")
parser.add_argument('-wavelet', action='store_true',  help="make wavelet plots for all variables")
parser.add_argument('-keogram', action='store_true', help="make keograms of all variables")
parser.add_argument('-heatmap', action='store_true', help="make 2D heat maps at different times, using r1,r2 as spatial variables")
parser.add_argument('-poynting', action='store_true', help="calculate the Poynting flux manually")
parser.add_argument('-cmap', default='plasma', help="colormap, for heatmaps and keograms")

args = parser.parse_args()

#run = 'EGP'    # 'EGL', 'EGP', ...  ALSO MODIFY virtual_sc.sh accordingly



fileprefix = args.fileprefix

if args.abs:
    fileprefix += '_abs'

if args.filter:
    fileprefix += '_filtered'

if args.dt:
    fileprefix += '_dt'

if args.norm:
    fileprefix += '_norm'



print('fileprefix ' + fileprefix)

#import pdb; pdb.set_trace()

root_dir = '/wrk-vakka/users/horakons/carrington/plots/'
save_dir = '{}{}/timeseries/'.format(root_dir, args.run.upper())         #CUSTOM


if args.n_sc is None:
    n_sc_input = None
else:
    n_sc_input= int(args.n_sc)


if args.filein is None:
    filein = 'virtual_sc_data_{}.txt'.format(args.run.lower())
else:
    filein = args.filein


print('reading file {}'.format(filein))
dct, nvars = read_vsc_data(args.filein, n_sc = n_sc_input, tvar = args.t)


#dct contains the original data as it was read from the file
#dct_plot will be cut down a bit (or even expanded upon), to use for analysis
dct_plot = deepcopy(dct)

# every item is a 2d array [i,j] where i indexes time (by default), j indexes position
# in the plots, however, t is the y-axis
t = dct_plot[args.t]
del dct_plot[args.t]

r1 = dct_plot[args.r1]
del dct_plot[args.r1]

try:
    r2 = dct_plot[args.r2]
    del dct_plot[args.r2]
except KeyError:
    r2 = r1*0
    print("r2 coordinate key error")
try:
    r3 = dct_plot[args.r3]
    del dct_plot[args.r3]
except KeyError:
    r3 = r1*0
    print("r3 coordinate key error")


global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations

global MU_0
MU_0 = 1.25663706e-6

global M_P
M_P = 1.67262192e-27

ntime = t.shape[0]
npos = t.shape[1]



if args.var is None:
    keys = sorted(dct_plot.keys())    # sort them alphabetically
else:
    # only consider the specified variables 
    keys = args.var

nvar = len(keys)






#filter data
#Savitzky-Golay Filter: Filter with a window length of X (must be odd) and a degree 2 polynomial. Use the defaults for all other parameters.
window_length = int(t.shape[0]/6)
if window_length % 2 == 0:
    window_length = window_length+1
#window_length = 201    # hard-code it (need to test that windowlength is smaller than array size, though)
for i in range(npos):
    print('location {}/{}'.format(i, npos))
    # make a time-series plot at the given coordinates of all the variables
    # test for vector components like B.x, B.y, B.z in a simple way: keep track of the last plotted thing,
    # if it was NAME.x and this is NAME.y (look for the '.') then don't close the current plot.
    # add a legend regardless of whether multiple things were plotted
    for j, key in enumerate(keys):
        try:
            if args.filter:
                filtered = savgol_filter(dct[key][:,i], window_length, 3)
                dct_plot[key][:,i] = dct[key][:,i]-filtered
        except KeyError:
            print("key {} not found. Skipping...".format(key))


#tack this on. (klug). S ~ dE x dB, which is why the filtering needs to happen first
if args.poynting:
    print(dct.keys())
    Sx = (1 / MU_0) * (dct_plot['vg_e_vol.y']*dct_plot['B.z'] - dct_plot['vg_e_vol.z']*dct_plot['B.y'])
    Sy = (1 / MU_0) * (dct_plot['vg_e_vol.z']*dct_plot['B.x'] - dct_plot['vg_e_vol.x']*dct_plot['B.z'])
    Sz = (1 / MU_0) * (dct_plot['vg_e_vol.x']*dct_plot['B.y'] - dct_plot['vg_e_vol.y']*dct_plot['B.x'])
    S = np.sqrt( Sx**2 + Sy**2 + Sz**2 )
    Spar = ( Sx*dct_plot['B.x'] + Sy*dct_plot['B.y'] + Sz*dct_plot['B.z'] ) / np.sqrt(dct_plot['B.x']**2 + dct_plot['B.y']**2 + dct_plot['B.z']**2)
    Sperp = np.sqrt( S**2 - Spar**2)
    dct_plot['Sx'] = Sx
    dct_plot['Sy'] = Sy
    dct_plot['Sz'] = Sz
    dct_plot['Spar'] = Spar
    dct_plot['Sperp'] = Sperp
    if args.var is None:
        keys = keys + ['Sx', 'Sy', 'Sz', 'Spar', 'Sperp']
        print('adding keys Sx, Sy, Sz, Spar, Sperp')
        keys = sorted(keys)    # maintain this list
        nvar = len(keys)


#more processing (normalize, d/dt, abs, etc.)
for i in range(npos):
    # make a time-series plot at the given coordinates of all the variables
    # test for vector components like B.x, B.y, B.z in a simple way: keep track of the last plotted thing,
    # if it was NAME.x and this is NAME.y (look for the '.') then don't close the current plot.
    # add a legend regardless of whether multiple things were plotted
    for j, key in enumerate(keys):
        try:
            if args.dt:
                dct_plot[0,i] = np.nan   #dummy 
                dct_plot[key][1:,i] = (dct_plot[key][1:,i] - dct_plot[key][0:ntime-1,i]) / (t[1:,i] - t[0:ntime-1,i])
            if args.norm:
                dct_plot[key][:,i] = dct_plot[key][:,i] / np.nanmax(np.abs(dct_plot[key][:,i]))
            if args.abs:
                dct_plot[key][:,i] = np.abs(dct_plot[key][:,i])
        except KeyError:
            print("key {} not found. Skipping...".format(key))





#now, make plots (timeseries, keograms, and/or heatmaps)
#timeseries
if args.timeseries:
    for i in range(npos):
        # make a time-series plot at the given coordinates of all the variables
        # add a legend regardless of whether multiple things were plotted
        fig, axes = plt.subplots(nvar, 1)      #axes is a 1D array with npos elements
        for j, key in enumerate(keys):
            label = make_label(key, filter=args.filter, dt=args.dt, abs = args.abs)
            axes[j].plot(t[:,i], dct_plot[key][:,i], label = label)
            if args.log:
                if np.min(dct_plot[key][:,i]) >= 0:
                    axes[j].set_yscale('log')
                else:
                    axes[j].set_yscale('symlog')
            axes[j].legend()
            if j < nvar-1:
                for xlabel_j in axes[j].get_xticklabels():
                    xlabel_j.set_visible(False)
            axes[j].set_ylim([str2float(args.vmin), str2float(args.vmax)])
        axes[j].xlabel = 'time [seconds]'
        axes[0].set_title( '{} virtual spacecraft ({},{},{})'.format(args.run.upper(), r1[0,i], r2[0,i], r3[0,i]) )
        plt.subplots_adjust(hspace=0)
        #fig.tight_layout()
        filename = '{}timeseries_multivar/{}_{}_coords_{}_{}_{}.png'.format(save_dir, fileprefix, args.run.upper(), r1[0,i], r2[0,i], r3[0,i])
        mkdir_path(filename)
        plt.savefig(filename)
        plt.close()



#wavelets
#pywavelets documentation: https://pywavelets.readthedocs.io/en/latest/ref/cwt.html
if args.wavelet:
    import pywt
    for i in range(npos):
        # make a wavelet plot at the given coordinates of all the variables
        # add a legend regardless of whether multiple things were plotted
        fig, axes = plt.subplots(nvar, 1)      #axes is a 1D array with npos elements
        widths = np.arange(1, 101)
        for j, key in enumerate(keys):
            label = make_label(key, filter=args.filter, dt=args.dt, abs = args.abs)
            cwtmatr, freqs = pywt.cwt(dct_plot[key][:,i], widths, 'mexh')
            filename = '{}timeseries_multivar/wavelet/{}_{}_{}_wavelet_coords_{}_{}_{}.png'.format(save_dir, fileprefix, args.run.upper(), key, r1[0,i], r2[0,i], r3[0,i])
            periods = 1/freqs
            #keogram.keogram(t[:,0], periods, cwtmatr, filename=filename, log = args.log, shading='auto', cmap = args.cmap,
            #                xlabel = 'time [sec], ' + ', ({},{},{}) to ({},{},{})',
            #                ylabel = 'frequency [Hz]', title = label + ', ' + args.run.upper() + ',  ({},{},{})'.format(r1[0,i], r2[0,i], r3[0,i]), cbar_label = label, 
            #                xlim = args.xlim, ylim = args.ylim, vmin = str2float(args.vmin), vmax = str2float(args.vmax))
            keogram.keogram(t[:,0], periods, cwtmatr, filename=filename, log = args.log, shading='auto', cmap = args.cmap,
                            xlabel = 'time [sec], ' + ', ({},{},{}) to ({},{},{})',
                            ylabel = 'period [sec]', title = label + ', ' + args.run.upper() + ',  ({},{},{})'.format(r1[0,i], r2[0,i], r3[0,i]), cbar_label = label, 
                            xlim = args.xlim, ylim = args.ylim, vmin = str2float(args.vmin), vmax = str2float(args.vmax))




#keograms
# e.g. cmap=plasma, vmin=0, vmax=1d9, shading='auto', etc.

#if args.keogram:
def keo(key):
    filename = '{}{}/keogram/{}/{}_keogram_{}_{}.png'.format(save_dir, key, fileprefix, fileprefix, key, args.run.upper())
    label = make_label(key, filter=args.filter, dt=args.dt, abs = args.abs)
    if args.dl:
        # add up the (approximate) length along the points. assume r1, r2, r3 are cartesian coordinates
        dx = r1[0,1:] - r1[0,0:npos-1]
        dy = r2[0,1:] - r2[0,0:npos-1]
        dz = r3[0,1:] - r3[0,0:npos-1]
        dl = (dx**2 + dy**2 + dz**2)**0.5
        l = np.zeros(npos)
        for i in range(npos-1):
            l[i+1] = l[i] + dl[i]
        keogram.keogram(l, t[:,0], dct_plot[key], filename=filename, log = args.log, shading='auto', cmap = args.cmap,
                        xlabel = 'L [R_E], ({},{},{}) to ({},{},{})'.format(r1[0,0], r2[0,0], r3[0,0], r1[0,npos-1], r2[0,npos-1], r3[0,npos-1]),
                        ylabel ='time [sec]', title = label + ', ' + args.run.upper(), cbar_label = label, vmin = str2float(args.vmin), vmax = str2float(args.vmax))   
    elif args.ta:
        # calculate distances in terms of the Alfven time dt_A = dl / v_A
        dx = r1[0,1:] - r1[0,0:npos-1]
        dy = r2[0,1:] - r2[0,0:npos-1]
        dz = r3[0,1:] - r3[0,0:npos-1]
        v_A = dct['vg_b_vol.magnitude'] / np.sqrt(MU_0 * dct['vg_rho'] * M_P)   # use dct instead of dct_plot, want to calculate v_A from unprocessed variables
        # compute the average alfven speed at each position
        v_A_ave = np.average(v_A, axis = 0)
        dl = (dx**2 + dy**2 + dz**2)**0.5
        t_A = np.zeros(npos)
        dt_A = 0
        for i in range(npos-1):
            if v_A_ave[i+1] != 0:
                dt_A = dl[i] * R_EARTH  / v_A_ave[i+1]
            t_A[i+1] = t_A[i] + dt_A
        keogram.keogram(t_A, t[:,0], dct_plot[key], filename=filename, log = args.log, shading='auto', cmap = args.cmap,
                        xlabel = 'Alfven time [sec], ' + ', ({},{},{}) to ({},{},{})'.format(r1[0,0], r2[0,0], r3[0,0], r1[0,npos-1], r2[0,npos-1], r3[0,npos-1]),
                        ylabel = 'time [sec]', title = label + ', ' + args.run.upper(), cbar_label = label, xlim = args.xlim, ylim = args.ylim, vmin = str2float(args.vmin), vmax = str2float(args.vmax))
    else:
        # every item in dct is a 2d array [x,y] where x indexes time, y indexes position
        keogram.keogram(r1[0,0:], t[:,0], dct_plot[key], filename=filename, log = args.log, shading='auto', cmap = args.cmap,
                        xlabel = args.r1.replace('_', '\\_') + ', ({},{},{}) to ({},{},{})'.format(r1[0,0], r2[0,0], r3[0,0], r1[0,npos-1], r2[0,npos-1], r3[0,npos-1]),
                        ylabel ='time [sec]', title = label + ', ' + args.run.upper(), cbar_label = label, xlim = args.xlim, ylim = args.ylim, vmin = str2float(args.vmin), vmax = str2float(args.vmax))


#heat maps (plot variables wrt r1, r2)
#for, e.g. making 2D visualizations of variables related to waves (ex. Spar)

r1_2d, r2_2d = squareAxes(r1[0,:], r2[0,:])

#if args.heatmap:
def heatmap(key):
    #first, detect if coordinates make up a square map
    vmin_temp = np.nanmin(dct_plot[key])
    vmax_temp = np.nanmax(dct_plot[key])
    vmax_abs = np.max([np.abs(vmin_temp), np.abs(vmax_temp)])
    if vmin_temp < 0:
        vmin = -vmax_abs
        vmax = vmax_abs
    else:
        vmin = vmin_temp
        vmax = vmax_temp
    # if set explicitly, override vmin, vmax
    if args.vmin is not None:
        vmin = float(args.vmin)
    if args.vmax is not None:
        vmax = float(args.vmax)
    print('vmin ' + str(vmin))
    print('vmax ' + str(vmax))
    label = make_label(key, filter=args.filter, dt=args.dt, abs=args.abs)
    for i in range(ntime):
        filename = '{}{}/heatmap/{}/{}_heatmap_{}_{}_{}.png'.format(save_dir, key, fileprefix, fileprefix, key, args.run.upper(), str(int(t[i,0])).zfill(5))
        print(filename)
        value = dct_plot[key][i,:].reshape(r1_2d.shape)
        keogram.keogram(r1_2d, r2_2d, value, filename=filename, log = args.log, shading='auto', cmap = args.cmap,
                        xlabel = args.r1.replace('_', '\\_'), ylabel = args.r2.replace('_', '\\_'), 
                        title = label + ', ' + args.run.upper() + ', t = ' + str(int(t[i,0])), cbar_label = label, xlim = args.xlim, ylim = args.ylim, vmin=vmin, vmax=vmax)


def parallelize_this(key):
    time_start = time()
    print(key)
    if args.keogram:
        keo(key)
    if args.heatmap:
        heatmap(key)
    time_stop = time()
    print('total time for key {} is: {} seconds'.format(key, time_stop - time_start))



## Parallel processing
from multiprocessing import Pool
pool = Pool(int(args.nproc))
return_array = pool.map(parallelize_this, keys)

pool.close()
pool.join()





