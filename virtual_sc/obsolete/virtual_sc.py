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



def poynting(Ex, Ey, Ez, Bx, By, Bz, filter=False, norm=False, window_length = 201):
    Ex_calc = deepcopy(Ex)
    Ey_calc = deepcopy(Ey)
    Ez_calc = deepcopy(Ez)
    Bx_calc = deepcopy(Bx)
    By_calc = deepcopy(By)
    Bz_calc = deepcopy(Bz)
    mu_0 = 1.25663706e-6
    for i in range(Ex.shape[1]):
        if filter:
            Ex_calc[:,i] = Ex[:,i] - savgol_filter(Ex[:,i], window_length, 3)
            Ey_calc[:,i] = Ey[:,i] - savgol_filter(Ey[:,i], window_length, 3)
            Ez_calc[:,i] = Ez[:,i] - savgol_filter(Ez[:,i], window_length, 3)
            Bx_calc[:,i] = Bx[:,i] - savgol_filter(Bx[:,i], window_length, 3)
            By_calc[:,i] = By[:,i] - savgol_filter(By[:,i], window_length, 3)
            Bz_calc[:,i] = Bz[:,i] - savgol_filter(Bz[:,i], window_length, 3)
    Sx = (1 / mu_0) * (Ey_calc*Bz_calc - Ez_calc*By_calc)
    Sy = (1 / mu_0) * (Ez_calc*Bx_calc - Ex_calc*Bz_calc)
    Sz = (1 / mu_0) * (Ex_calc*By_calc - Ey_calc*Bx_calc)
    S = np.sqrt( Sx**2 + Sy**2 + Sz**2 )
    Spar = ( Sx*Bx + Sy*By + Sz*Bz ) / np.sqrt(Bx**2 + By**2 + Bz**2)
    Sperp = np.sqrt( S**2 - Spar**2)
    if norm:
        Sx = Sx / S
        Sy = Sy / S
        Sz = Sz / S
        Spar = Spar / S
        Sperp = Sperp / S
    return Sx, Sy, Sz, Spar, Sperp


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
parser.add_argument('-xlim', nargs='*', help="axes xlim (list [xmin, xmax])" )
parser.add_argument('-ylim', nargs='*', help="axes ylim (list [xmin, xmax])" )
#parser.add_argument('-vmin', nargs='*', help="list of mins (one for each variable)" )
#parser.add_argument('-vmax', nargs='*', help="list of maxes (one for each variable)" )
parser.add_argument('-filter', action='store_true',  help="apply Savitsky-Golay filter to all data, subtract filtered data to get delta")
parser.add_argument('-log', action='store_true',  help="plot logarithm of z-data")
parser.add_argument('-norm', action='store_true',  help="if filter is set, divide filtered data by trend. if filter is not set, divide by background (initial condition)")
parser.add_argument('-timeseries', action='store_true',  help="make time series plots of all variables")
parser.add_argument('-keogram', action='store_true', help="make keograms of all variables")
parser.add_argument('-heatmap', action='store_true', help="make 2D heat maps at different times, using r1,r2 as spatial variables")
parser.add_argument('-poynting', action='store_true', help="calculate the Poynting flux manually")
args = parser.parse_args()

#run = 'EGP'    # 'EGL', 'EGP', ...  ALSO MODIFY virtual_sc.sh accordingly



fileprefix = args.fileprefix

if args.filter:
    fileprefix += '_filtered'

if args.norm:
    fileprefix += '_norm'


print('fileprefix ' + fileprefix)

#enable debugging
#import pdb; pdb.set_trace()

root_dir = '/wrk-vakka/users/horakons/carrington/plots/'
save_dir = '{}{}/timeseries/'.format(root_dir, args.run.upper())         #CUSTOM


if args.n_sc is None:
    n_sc_input = None
else:
    n_sc_input= int(args.n_sc)


if args.filein is None:
    dct, nvars = read_vsc_data('txt_files/virtual_sc_data_{}.txt'.format(args.run.lower()), n_sc = n_sc_input, tvar=args.t)
else:
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


ntime = t.shape[0]
npos = t.shape[1]



if args.var is None:
    keys = sorted(dct_plot.keys())    # sort them alphabetically
else:
    # only consider the specified variables 
    keys = args.var

nvar = len(keys)




#Savitzky-Golay Filter: Filter with a window length of X (must be odd) and a degree 2 polynomial. Use the defaults for all other parameters.
window_length = int(t.shape[0]/6)
if window_length % 2 == 0:
    window_length = window_length+1
#window_length = 201    # hard-code it (need to test that windowlength is smaller than array size, though)


#process data (filter, normalize, etc.)
for i in range(npos):
    print('location {}/{}'.format(i, npos))
    # make a time-series plot at the given coordinates of all the variables
    # test for vector components like B.x, B.y, B.z in a simple way: keep track of the last plotted thing,
    # if it was NAME.x and this is NAME.y (look for the '.') then don't close the current plot.
    # add a legend regardless of whether multiple things were plotted
    for j, key in enumerate(keys):
        try:
            value = dct[key]
            if args.filter:
                filtered = savgol_filter(value[:,i], window_length, 3)
                if args.norm:
                    dct_plot[key][:,i] = (value[:,i]-filtered) / filtered
                else:
                    dct_plot[key][:,i] = value[:,i]-filtered
            else:
                if args.norm:
                    dct_plot[key][:,i] = value[:,i] / np.abs(value[0,i])
                else:
                    dct_plot[key][:,i] = value[:,i]
        except KeyError:
            print("key {} not found. Skipping...".format(key))


#tack this on. (klug). Need to process the poynting flux manually
if args.poynting:
    print(dct.keys())
    Sx, Sy, Sz, Spar, Sperp = poynting( dct['vg_e_vol.x'], dct['vg_e_vol.y'], dct['vg_e_vol.z'], dct['B.x'], dct['B.y'], dct['B.z'],
                                        filter = args.filter, norm=args.norm, window_length = window_length )
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


#now, make plots (timeseries, keograms, and/or heatmaps)
#timeseries
if args.timeseries:
    for i in range(npos):
        # make a time-series plot at the given coordinates of all the variables
        # add a legend regardless of whether multiple things were plotted
        fig, axes = plt.subplots(nvar, 1)      #axes is a 1D array with npos elements
        for j, key in enumerate(keys):
            if args.filter:
                label = r"$\Delta$ " + key.replace('_', '\\_')
            else:
                label = key.replace('_', '\\_')
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
        axes[j].xlabel = 'time [seconds]'
        plt.subplots_adjust(hspace=0)
        #fig.tight_layout()
        filename = '{}{}_{}_coords_{}_{}_{}.png'.format(save_dir, fileprefix, args.run.upper(), r1[0,i], r2[0,i], r3[0,i])
        mkdir_path(filename)
        plt.savefig(filename)
        plt.close()



#keograms
# e.g. cmap=plasma, vmin=0, vmax=1d9, shading='auto', etc.

#if args.keogram:
def keo(key):
    filename = '{}{}_keogram_{}_{}.png'.format(save_dir, fileprefix, key, args.run.upper())
    value = dct_plot[key]
    if args.filter:
        label = r"$\Delta$ " + key.replace('_', '\\_')
    else:
        label = key.replace('_', '\\_')
    if args.dl:
        # add up the (approximate) length along the points. assume r1, r2, r3 are cartesian coordinates
        dx = r1[0,1:] - r1[0,0:npos-1]
        dy = r2[0,1:] - r2[0,0:npos-1]
        dz = r3[0,1:] - r3[0,0:npos-1]
        dl = (dx**2 + dy**2 + dz**2)**0.5
        l = np.zeros(npos)
        for i in range(npos-1):
            l[i+1] = l[i] + dl[i]
        keogram.keogram(l, t[:,0], value, filename=filename, log = args.log, shading='auto',
                        xlabel = 'L [AU], ({},{},{}) to ({},{},{})'.format(r1[0,0], r2[0,0], r3[0,0], r1[0,npos-1], r2[0,npos-1], r3[0,npos-1]),
                        ylabel ='time [sec]', title = label + ', ' + args.run.upper(), cbar_label = label)   
    else:
        # every item in dct is a 2d array [x,y] where x indexes time, y indexes position
        keogram.keogram(r1[0,0:], t[:,0], value, filename=filename, log = args.log, shading='auto',
                        xlabel = args.r1.replace('_', '\\_'),
                        ylabel ='time [sec]', title = label + ', ' + args.run.upper(), cbar_label = label, xlim = args.xlim, ylim = args.ylim)


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
    if args.filter:
        label = r"$\Delta$ " + key.replace('_', '\\_')
    else:
        label = key.replace('_', '\\_')
    for i in range(ntime):
        filename = '{}{}/{}_heatmap_{}_{}_{}.png'.format(save_dir, key, fileprefix, key, args.run.upper(), str(int(t[i,0])).zfill(5))
        value = dct_plot[key][i,:].reshape(r1_2d.shape)
        if args.filter:
            label = r"$\Delta$ " + key.replace('_', '\\_')
        else:
            label = key.replace('_', '\\_')
        keogram.keogram(r1_2d, r2_2d, value, filename=filename, log = args.log, shading='auto',
                        xlabel = args.r1.replace('_', '\\_'), ylabel = args.r2.replace('_', '\\_'), 
                        title = label + ', ' + args.run.upper() + ', t = ' + str(int(t[i,0])), cbar_label = label, xlim = args.xlim, ylim = args.ylim, vmin=vmin, vmax=vmax)


if __name__=='__main__':
    for key in keys:
        time_start = time()
        print(key)
        if args.keogram:
            keo(key)
        if args.heatmap:
            heatmap(key)
        time_stop = time()
        print('total time for key {} is: {} seconds'.format(key, time_stop - time_start))



### Parallel processing
#from multiprocessing import Pool
#pool = Pool(int(args.nproc))
#return_array = pool.map(parallelize_this, keys)






