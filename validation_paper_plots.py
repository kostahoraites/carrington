# make plots for the Vlasiator pressure pulse validation paper

import pytools as pt
from pyCalculations.intpol_points import vlsv_intpol_points
import pyCalculations.fieldtracer as fieldtracer
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from static_field_tracer_3d_alt import static_field_tracer_3d_alt
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path, cartesian_to_spherical, spherical_to_cartesian, numcurl3d, numjacobian3d
import os, sys
import argparse
import warnings
from copy import deepcopy
import scipy
import pandas as pd
from time import time
from carrington import *
from tsyganenko import *
from cusp import cusp

mpl.rcParams.update({'font.size': 18})   # increase font size


# example call from bash (see carrington.sh):
# python carrington.py -run EGI -var 1 4 6

# Input parameters
parser = argparse.ArgumentParser()

parser.add_argument('-nproc', default=1, help="number of processors to use " )
global ARGS
ARGS = parser.parse_args()
#nproc = ARGS.nproc


global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations
global plot_dir
global mu_0
mu_0 = 4e-7 * np.pi

#run = 'EGL'   # EGP
plot_dir = '/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/'

########### Figure 1a) 3D perspective plot of density ##########



#EXAMPLE 5: plot_threeslice

import pytools as pt


########## Figure 1b) Pearson correlation, best fit with OMNI ###########


########## Figure 2 (setup) ###########



def f_ocb_control(phi_i):
    return tsyganenko_ocb(phi_i, lat_range = [60,80], nsteps = 10,
                          Txx = 't01', Dst = -33, Kp = 4, Vx_sw = -750., N_sw = 1., Bx_imf = 0., By_imf = 0., Bz_imf = -5.,
                          R_inner = 1., R_outer=15., dir = None, maxloop = 10000)


def f_ocb_pulse(phi_i):
    return tsyganenko_ocb(phi_i, lat_range = [60,80], nsteps = 10,
                          Txx = 't01', Dst = -33, Kp = 4, Vx_sw = -750., N_sw = 4., Bx_imf = 0., By_imf = 0., Bz_imf = -10.,
                          R_inner = 1., R_outer=15., dir = None, maxloop = 10000)


# define quantities normally contained in ARGS
    
phimin = -180.
phimax = 180.
latmin = 60.  # (CHANGE THIS!)
latmax = 90.
nlat = 301   #121
nphi = 361   #361  

# make longitude-latitude grid over which to evaluate data
#nlat = 180
#nphi = 360
#lat, phi = lat_phi_grid( nlat = nlat, nphi = nphi)

degtorad = np.pi / 180
lat, phi = lat_phi_grid( phi_min = float(phimin)*degtorad, phi_max = float(phimax)*degtorad,
                         lat_min = float(latmin)*degtorad, lat_max = float(latmax)*degtorad,
                         nlat = int(nlat), nphi = int(nphi))


#calculate tsyganenko OCB boundary

phi_azim_plot = phi[0,:] + (np.pi / 2)   # rotate by 90 degrees so noon is at the top
phi_1d = phi[0,:] * 180 / np.pi
lat_1d = lat[:,0] * 180 / np.pi

## v   Parallel processing     v

#UNCOMMENT TO RECALCULATE TSYGANENKO OCB

#from multiprocessing import Pool
#pool = Pool(int(ARGS.nproc))
#ocb_control = pool.map(f_ocb_control, phi_1d)
#pool.close()
#pool.join()

#pool = Pool(int(ARGS.nproc))
#ocb_pulse = pool.map(f_ocb_pulse, phi_1d)
#pool.close()
#pool.join()

#np.savetxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_control.txt', ocb_control)    
#np.savetxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_pulse.txt', ocb_pulse)

## ^        ^

ocb_tsyg_control = np.loadtxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_before.txt')    
ocb_tsyg_pulse = np.loadtxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_after.txt')
phi_1d_tsyg = np.linspace(-np.pi, np.pi, num = ocb_tsyg_pulse.size)
phi_azim_tsyg_plot = phi_1d_tsyg + (np.pi / 2)   # rotate by 90 degrees so noon is at the top

phi_mlt = phi_1d * (24. / 360) + 12
phi_mlt_tsyg = phi_1d_tsyg * (12. / np.pi) + 12



##

def ocb_2d_to_1d(ocb_2d,lat_1d):
    nlat = ocb_2d.shape[0]
    nphi = ocb_2d.shape[1]
    ocb_1d = np.zeros(nphi)
    for i in range(nphi):
        latcut = ocb_2d[:,i]
        clsd_init = latcut[0]
        #for j in range(1, nlat):               #find the most poleward topology change
        for j in reversed(range(1, nlat)):      #find the equatormost topology change
            if latcut[j] != latcut[j-1]:  # find topology changes, closed <--> open 
                 ocb_1d[i] = ( (lat_1d[j] + lat_1d[j-1]) / 2.)  # already in degrees
    return ocb_1d





data = {'Control':{'fileIndex':-1, 'ocb_tsyg':ocb_tsyg_control}, 'Pulse':{'fileIndex':-1, 'ocb_tsyg':ocb_tsyg_pulse}}
#data = {'Control':{'fileIndex':-1, 'ocb_tsyg':ocb_tsyg_control}}

# aurora plots:

labels = ['a)', 'c)', 'e)', 'b)', 'd)', 'f)']
ilabel = 0


for key in data: # 'Control', 'Pulse'

    # Defining source and output file locations

    if key == 'Control':
        run = 'EGI'   # EGP
        j=857                       #(j here is from EGL arrival time)
        # EGI: 662-1506
        data[key]['fileIndex'] = 1506
    if key == 'Pulse':
        run = 'EGL'
        j=857                       #EGL --- pressure pulse arrival time
        #     EGL: 621 - 1760
        data[key]['fileIndex'] = 1760
    #dim = '3D'
    #bulk = 'bulk1.egl'
    #bulkLocation = '/wrk/group/spacephysics/vlasiator/{}/{}/bulk/'.format(dim, run)
    #fluxLocation = '/wrk/group/spacephysics/vlasiator/{}/{}/flux/'.format(dim, run)
    outputLocation = '/wrk-vakka/users/horakons/carrington/plots/{}_validation_paper/'.format(run)
    #bulkname = '{}.{}.vlsv'.format(bulk,str(j).zfill(7))
    #elif run == 'EGP':
    #    j=506                       #EGP --- first time with precipitation data   (269-506 total)
    #    control_index = 352
    #    pulse_index = 506

    bulkname = get_vlsvfile_fullpath(run, j)

    ##pt.plot.plot_threeslice(filename=bulkLocation+bulkname, var ='proton/vg_rho', colormap = 'plasma', vmin = 1e4,vmax=5e6, step=j, outputdir=outputLocation,outputfile='test_threeslice_{}_{}.jpg'.format(run,j), Earth=1,cutpointre=[0,0,0], slices='yz' )

    '''
    pt.plot.plot_threeslice(filename=bulkname, var ='vg_pdyn', colormap = 'plasma', vmin = 1e-10,vmax=5e-9, step=j, outputdir=outputLocation,outputfile='test_threeslice_pdyn_{}_{}.jpg'.format(run,j), Earth=1,cutpointre=[0,0,0], slices='yz', scale=1e-9, usesci=0 )
    '''

    #plt.close()

    # load data
    fileIndex = data[key]['fileIndex']
    filename = get_vlsvfile_fullpath(run, fileIndex)
    f = pt.vlsvfile.VlsvReader( filename )
    r_aurora = 1. * R_EARTH + 1.1e5    # auroral altitude 110 km
    r_trace = 5 * R_EARTH
    dx = R_EARTH / 50
    max_iterations = int(4 * r_trace / dx)

    fg_b = f.read_variable('fg_b')     # this is SLOW  ... TODO: implement hongyang's julia version  

    data[key]['fg_b'] = fg_b

    global CELLSIZE_XYZ
    CELLSIZE_XYZ = [ (f.read_parameter('xmax') - f.read_parameter('xmin')) / fg_b.shape[0],
                    (f.read_parameter('ymax') - f.read_parameter('ymin')) / fg_b.shape[1],
                    (f.read_parameter('zmax') - f.read_parameter('zmin')) / fg_b.shape[2] ]
    
    theta = lat2theta(lat)
    coord_list_aurora = [(phi*0) + r_aurora, theta, phi]
    coord_list_aurora_cart = spherical_to_cartesian( *coord_list_aurora )
    coord_list_magnetosphere = [(phi*0) + r_trace, theta, phi]
    coord_list_magnetosphere_cart = spherical_to_cartesian( *coord_list_magnetosphere )

    # trace lines in a radial orientation
    
    #trace_coord = sidecar(trace_coord)        # SIDECAR: read function output from hard disk if calulated control
    
    trace_method = 'integrateB'            #options: 'dipole', 'integrateB'
    coord_list_traced_pos = trace_coord(filename, coord_list_aurora, fg_b = fg_b, trace_method = trace_method, direction = '+',
                                        coord_in = 'spherical', coord_out = 'spherical', r_trace = r_trace, max_iterations = max_iterations, dx = dx)
    coord_list_traced_neg = trace_coord(filename, coord_list_aurora, fg_b = fg_b, trace_method = trace_method, direction = '-',
                                        coord_in = 'spherical', coord_out = 'spherical', r_trace = r_trace, max_iterations = max_iterations, dx = dx)
    mask = test_radial_field(f, coord_list_aurora_cart)
    coord_list_traced = coord_list_traced_neg
    for i in range(3):
        coord_list_traced[i][mask] = coord_list_traced_pos[i][mask]
    
    # make plots and save the data for different parameters in objects 
    ########### Figure 2a) open/closed field lines ##########
    # open vs. closed field boundary

    tsyg_label = 'T01'
    ocb_label = 'OCB'

    fOO = FieldOpenObj(run, f, coord_list_aurora, fg_b = fg_b, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data = False)
    ax = fOO.make_plot_data( plot_data=True, plots = [0], pltclose=False)  # polar plot
    #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
    ax.set_title('OCB, ({}), t={} sec'.format(key, fileIndex), pad = 16)
    plt.text(0.05, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 24)
    ilabel += 1
    #plt.legend()
    plt.subplots_adjust(top=0.85)
    plt.savefig(plot_dir+'ocb_{}_{}_azim.pdf'.format(key, run))
    plt.close()

    data[key]['ocb_2d_tf'] = fOO.data_list[0]     # boolean, [nlat, nphi]

    data[key]['ocb_1d'] = ocb_2d_to_1d(data[key]['ocb_2d_tf'], lat_1d)
    data[key]['ocb_1d'] = ocb_2d_to_1d(data[key]['ocb_2d_tf'], lat_1d)

    ########### Figure 2b) FACs ##########

    # Field aligned currents (FACS)

    JO_traced = FAC_Obj(run, f, coord_list_aurora, fg_b = fg_b, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data = False)
    ax = JO_traced.make_plot_data( plot_data=True, plots = [1], pltclose=False)  # polar plot 
    plt.sca(ax) 
    #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
    ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label=ocb_label, color = "green" )
    plt.text(0.05, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 24)
    ilabel += 1
    ax.set_title('FAC, ({}), t={} sec'.format(key, fileIndex), pad = 16)
    #plt.legend()

    
    # overplot MFACE model

    '''
    mface_dir = '/wrk-vakka/users/horakons/carrington/MFACE/'
    if key == 'Control':
        outpt = pd.read_csv(mface_dir+'output_before.txt', delim_whitespace = True)
    elif key == 'Pulse':
        outpt = pd.read_csv(mface_dir+'output_after.txt', delim_whitespace = True)
    inpt = pd.read_csv(mface_dir+'input.txt', delim_whitespace = True)

    MLT = np.array(inpt['MLT(hour)'])
    MLAT = np.array(inpt['MLAT(degree)'])
    J = np.array(outpt['MeanJ(uA/m2)'])


    #fig = plt.gcf()
    #gs = gridspec.GridSpec(1, 2, width_ratios=[10,1])
    ##gs = gridspec.GridSpec(2, 1)
    #ax2 = plt.subplot(gs[0], projection="polar", aspect=1.)
    ##phi_plot = phi + (np.pi / 3)          # rotate by ninety-degrees so noon is at the top
    ##theta_plot = theta*181/np.pi

    MLT = np.array(inpt['MLT(hour)'])
    MLAT = np.array(inpt['MLAT(degree)'])
    J = np.array(outpt['J(uA/m2)'])

    phi_plot = MLT * 361 / 24.
    theta_plot = 91. - MLAT 
    phi_plot = phi_plot.reshape([360,121]) - 90. 
    theta_plot = theta_plot.reshape([360,121])
    J = J.reshape([360,121])

    if np.nanmin(theta_plot) > 91:
        theta_plot = 181 - theta_plot
    #im = ax2.pcolormesh(phi_plot*np.pi / 180, np.abs(theta_plot), J,
    #                    cmap='bwr', shading='auto')            # north pole
    print('plotting contour...')
    im = ax.contour(phi_plot*np.pi / 180, np.abs(theta_plot), J, 7, cmap='bwr', vmin = -2., vmax = 2.);
    #im = ax2.contour(phi_plot*np.pi / 180, np.abs(theta_plot), J, 13, cmap='PuOr', vmin = -3., vmax = 3.);
    #ax2.plot(phi_plot.flatten(), phi_plot.flatten() * 0 + 10)
    #ax2.set_ylim([0,40])
    #ax2.set_yticks([10, 20, 30])  # Less radial ticks
    #ax2.set_yticklabels([r'$80^\circ$', r'$70^\circ$', r'$60^\circ$'])
    #ax2.set_rlabel_position(-10)  # Move radial labels away from plotted line
    #ax2.set_xticks(list(np.arange(0, 2*np.pi, np.pi / 4)))
    #ax.set_ylim([0,40])
    #ax.set_yticks([10, 20, 30])  # Less radial ticks
    #ax.set_yticklabels([r'$80^\circ$', r'$70^\circ$', r'$60^\circ$'])
    #ax.set_rlabel_position(-10)  # Move radial labels away from plotted line
    #ax.set_xticks(list(np.arange(0, 2*np.pi, np.pi / 4)))
    #ax2.set_xticklabels(list(np.arange(6, 24, 3)) + list(np.arange(0, 6, 3)))
    #ax2.grid(True)
    #ax.set_xticklabels(list(np.arange(6, 24, 3)) + list(np.arange(0, 6, 3)))
    #ax.grid(True)
    #ax3 = plt.subplot(gs[1])
       #cax = fig.add_axes([ax.get_position().x2+0.01,ax.get_position().y0,0.02,ax.get_position().height])
       #plt.colorbar(im, cax=cax) # Similar to fig.colorbar(im, cax = cax)
       #divider = make_axes_locatable(ax2)
       #cax = divider.append_axes('right', size='6%', pad=0.05)
       #cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    #cbar = plt.colorbar(im, cax=ax3, aspect = 40)

    #divider = make_axes_locatable(ax)                  # tried this, didn't work
    #cax = divider.append_axes('left', size='5%', pad=0.05)
    #mycbar = plt.colorbar(im, cax=cax, orientation='vertical')
    #mycbar.set_label(r'J [$\mu A/m^3$]')

    #mycbar = fig.colorbar(im, cax=cax, orientation='vertical')
    #mycbar = ax.figure.colorbar(im)

    #cbar = plt.colorbar(im, cax=ax, aspect = 40)
       #cbar = mpl.colorbar.ColorbarBase(ax3, cmap = cmap)
    #cbar.set_label(r'J [$\mu A/m^3$]')
    #ax2.set_title('MFACE FAC')
    '''

    # OVERPLOT THE pyAMPS model
    
    mface_dir = '/wrk-vakka/users/horakons/carrington/pyAMPS/'
    if key == 'Control':
        outpt = pd.read_csv(mface_dir+'pyamp_Control.csv', delim_whitespace = True)
    elif key == 'Pulse':
        outpt = pd.read_csv(mface_dir+'pyamp_Pulse.csv', delim_whitespace = True)

    mask = (np.array(outpt['mlat']) > 0)      # 10,000 elements in northern hemisphere, mlat> 0

    MLT = np.array(outpt['mlt'])[mask]
    MLAT = np.array(outpt['mlat'])[mask]
    J = np.array(outpt['Jd'])[mask]    # units?

    fig = plt.gcf()
    #gs = gridspec.GridSpec(1, 2, width_ratios=[10,1])
    ##gs = gridspec.GridSpec(2, 1)
    #ax2 = plt.subplot(gs[0], projection="polar", aspect=1.)
    ##phi_plot = phi + (np.pi / 3)          # rotate by ninety-degrees so noon is at the top
    ##theta_plot = theta*181/np.pi

    phi_plot = MLT * 361 / 24.
    theta_plot = 91. - MLAT 
    phi_plot = phi_plot.reshape([100,100]) - 90. 
    theta_plot = theta_plot.reshape([100,100])
    J = J.reshape([100,100])

    if np.nanmin(theta_plot) > 91:
        theta_plot = 181 - theta_plot
    #im = ax2.pcolormesh(phi_plot*np.pi / 180, np.abs(theta_plot), J,
    #                    cmap='bwr', shading='auto')            # north pole
    print('plotting contour...')
    #im = ax.contour(phi_plot*np.pi / 180, np.abs(theta_plot), J, [-2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5], cmap='PuOr', vmin = -2.5, vmax = 2.5)   # linestyles= 'dashed'
    im = ax.contour(phi_plot*np.pi / 180, np.abs(theta_plot), J, [-2.4, -1.6, -0.8, 0.8, 1.6, 2.4], cmap='PuOr', vmin = -2.5, vmax = 2.5)   # linestyles= 'dashed'
    #im = ax2.contour(phi_plot*np.pi / 180, np.abs(theta_plot), J, 13, cmap='PuOr', vmin = -3., vmax = 3.);
    #ax2.plot(phi_plot.flatten(), phi_plot.flatten() * 0 + 10)
    #ax2.set_ylim([0,40])
    #ax2.set_yticks([10, 20, 30])  # Less radial ticks
    #ax2.set_yticklabels([r'$80^\circ$', r'$70^\circ$', r'$60^\circ$'])
    #ax2.set_rlabel_position(-10)  # Move radial labels away from plotted line
    #ax2.set_xticks(list(np.arange(0, 2*np.pi, np.pi / 4)))
    #ax.set_ylim([0,40])
    #ax.set_yticks([10, 20, 30])  # Less radial ticks
    #ax.set_yticklabels([r'$80^\circ$', r'$70^\circ$', r'$60^\circ$'])
    #ax.set_rlabel_position(-10)  # Move radial labels away from plotted line
    #ax.set_xticks(list(np.arange(0, 2*np.pi, np.pi / 4)))
    #ax2.set_xticklabels(list(np.arange(6, 24, 3)) + list(np.arange(0, 6, 3)))
    #ax2.grid(True)
    #ax.set_xticklabels(list(np.arange(6, 24, 3)) + list(np.arange(0, 6, 3)))
    #ax.grid(True)
    #ax3 = plt.subplot(gs[1])
       #cax = fig.add_axes([ax.get_position().x2+0.01,ax.get_position().y0,0.02,ax.get_position().height])
       #plt.colorbar(im, cax=cax) # Similar to fig.colorbar(im, cax = cax)
       #divider = make_axes_locatable(ax2)
       #cax = divider.append_axes('right', size='6%', pad=0.05)
       #cbar = fig.colorbar(im, cax=cax, orientation='vertical')
    #cbar = plt.colorbar(im, cax=ax3, aspect = 40)

    #divider = make_axes_locatable(ax)                  # tried this, didn't work
    #cax = divider.append_axes('left', size='5%', pad=0.05)
    #mycbar = plt.colorbar(im, cax=cax, orientation='vertical')
    #mycbar.set_label(r'J [$\mu A/m^3$]')

    class nf(float):
        def __repr__(self):
            s = f'{self:.1f}'
            return f'{self:.0f}' if s[-1] == '0' else s

    # Recast levels to new class
    im.levels = [nf(val) for val in im.levels]

    # Label levels with specially formatted floats
    if plt.rcParams["text.usetex"]:
        fmt = r'%r'
    else:
        fmt = '%r'

    ax.clabel(im, im.levels, inline=True, fmt=fmt, fontsize=10, colors = 'black')

    #fig.colorbar(im, cax=cax, orientation='vertical')
    #mycbar = ax.figure.colorbar(im)

    #cbar = plt.colorbar(im, cax=ax, aspect = 40)
       #cbar = mpl.colorbar.ColorbarBase(ax3, cmap = cmap)
    #cbar.set_label(r'J [$\mu A/m^3$]')
    #ax2.set_title('MFACE FAC')
    

    plt.subplots_adjust(top=0.85)
    print(plot_dir+'fac_{}_{}_azim.pdf'.format(key, run))
    plt.savefig(plot_dir+'fac_{}_{}_azim.pdf'.format(key, run))
    plt.close()

    ########### Figure 2c) precipitation ##########

    # proton differential energy flux
    #if run == 'EGL':   # EGP doesn't have precipitation data control t=352
    if True:   # EGP doesn't have precipitation data
        i0 = 5
        #if run == 'EGL':
        #    i0 = 5
        #elif run == 'EGP':
        #    i0 = 9
        pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data=False)
        i = i0
        ax = pFO_traced.make_plot_data( plot_data=True, plots = [i], pltclose=False)  # polar plot
        #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
        ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label = ocb_label, color = "green" )
        plt.text(0.05, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 24)
        proton_energy = f.read_parameter('proton_PrecipitationCentreEnergy{}'.format(i))
        ax.set_title('DNF, {:.1f} keV, ({}), t={} sec'.format(proton_energy/1000., key, fileIndex), pad = 16)
        #plt.legend()  '{:.1f}'.format(8891 / 1000.)
        plt.subplots_adjust(top=0.85)
        plt.savefig(plot_dir+'dnf_{}_{}_{}_azim.pdf'.format(key, int(proton_energy), run))
        plt.close()
    
        pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data=False)
        i = i0+1
        ax = pFO_traced.make_plot_data( plot_data=True, plots = [i], pltclose=False)  # polar plot
        #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
        ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label = ocb_label, color = "green" )
        plt.text(0.05, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 24)
        ilabel += 1
        proton_energy = f.read_parameter('proton_PrecipitationCentreEnergy{}'.format(i))
        ax.set_title('DNF, {:.1f} keV, t={} sec'.format(proton_energy/1000., fileIndex), pad = 16)
        #plt.legend()
        plt.subplots_adjust(top=0.85)
        plt.savefig(plot_dir+'dnf_{}_{}_{}_azim.pdf'.format(key, int(proton_energy), run))
        plt.close()
    

        # proton differential energy flux

        #pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data=False)
        #ax = pFO_traced.make_plot_data( plot_data=True, plots = [5], pltclose=False)  # polar plot
        ##ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
        #ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label = ocb_label, color = "green" )
        #plt.text(0.1, 0.9, labels[ilabel], transform=plt.gcf().transFigure)
        #ilabel += 1
        #ax.set_title('DNF, t={} sec'.format(fileIndex))
        ##plt.legend()
        #plt.savefig(plot_dir+'def_{}_azim.pdf'.format(key))
        #plt.close()

        #pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = True, plot_data=plot)
        #save_list.append(pFO_traced)



# Figure 4: OCB vs. MLT   (with Tsyganenko?)

ocb_1d_control = data['Control']['ocb_1d']
ocb_1d_pulse = data['Pulse']['ocb_1d']
#ocb_1d_control = data['control']['ocb_tsyg']
#ocb_1d_pulse = data['pulse']['ocb_tsyg']

#plt.plot( phi_azim_plot, ocb_1d_control)
#plt.savefig(plot_dir+'ocb_control_line.png')
#plt.close()

#plt.plot( phi_azim_plot, ocb_1d_pulse)
#plt.savefig(plot_dir+'ocb_pulse_line.png')
#plt.close()

fig, ax = plt.subplots()
plt.xlabel(r'MLT [hours]')
plt.ylabel(r'OCB Latitude [$^{\circ}$]')
plt.title(r'OCB vs. MLT')
plt.plot(phi_mlt, ocb_1d_control, label='OCB (Control)' )
plt.plot(phi_mlt, ocb_1d_pulse, label='OCB (Pulse)' )
plt.ylim([60, 79])
#plt.ylim([72, 78])
ax.set_xticks([6,8,10,12,14,16,18])
ax.set_xlim([6,18])
#ax.set_xticks([0,4,8,12,16,20,24])
#plt.plot(phi_mlt_tsyg, ocb_tsyg_control, label='Tsyg. (control)' )
#plt.plot(phi_mlt_tsyg, ocb_tsyg_pulse, label='Tsyg. (pulse)' )
#plt.xlim([0,24])

# Compare with newell et al., 2006 (Table 1) for the cusp latitude

names =  ['Bz', 'EKL', 'Epsilon', 'Rectifier', 'Pressure', 'Vasyliunas', 'Density', 'Velocity', 'WAV']
plot_names =  [r'$B_z$', r'$E_{KL}$', r'$\epsilon$', 'Rectifier', 'pressure', 'Vasyliunas', 'density', 'velocity', 'WAV']

# 'Control'
v = 750.
nn = 1.
p = 1.67e-27 * (nn * 1e6) * (v * 1e3) * 1e9
B = 5.
Bz = -5
Bt = 5.
Bs = np.min([0., Bz])
theta_c = np.pi
for n, pn in zip(names, plot_names):
    c = cusp(n, v = v, n = nn, p = p, B = B, Bz = Bz, Bt = Bt, Bs = Bs, theta_c = theta_c)
    print(n, c)
    ax.scatter(11.9, c, color = 'tab:blue', facecolors = 'none')
    if n=='Vasyliunas':
        xoffset = 0
        yoffset = -0.4
    elif n=='WAV':
        xoffset = 0
        yoffset = 0.4
    else:
        xoffset=0; yoffset =0
    ax.annotate(pn, (11.5+xoffset, c+yoffset), color = 'tab:blue', fontsize = 10, horizontalalignment = 'right')

# 'Pulse'
v = 750.
nn = 4.
p = 1.67e-27 * (nn * 1e6) * (v * 1e3) * 1e9
B = 10.
Bz = -10
Bt = 10.
Bs = np.min([0., Bz])
theta_c = np.pi
for n, pn in zip(names, plot_names):
    c = cusp(n, v = v, n = nn, p = p, B = B, Bz = Bz, Bt = Bt, Bs = Bs, theta_c = theta_c)
    print(n, c)
    ax.scatter(12.1, c, color = 'tab:orange', facecolors = 'none')
    if n=='Vasyliunas':
        xoffset = 0
        yoffset = -0.4
    elif n=='WAV':
        xoffset = 0
        yoffset = -0.4
    else:
        xoffset=0; yoffset =0
    ax.annotate(pn, (12.5+xoffset, c+yoffset), color = 'tab:orange', fontsize = 10, horizontalalignment = 'left')

plt.legend(loc = 'lower left', fontsize = 'small')
plt.tight_layout() 
plt.savefig(plot_dir + 'ocb_vs_mlt.pdf')
plt.savefig(plot_dir + 'ocb_vs_mlt.png')
plt.close()

#RESULT: data has 'fileIndex', 'ocb_2d_tf', and 'fg_b' keys






########### Figure 3a) Magnetopause position R from beta*, x-y plane ##########
########### Figure 3b) Magnetopause position R from beta*, x-z plane ##########

# see magnetopause/run_beta_star_colormap.sh   # fiddle with the time step as necessary


########### Figure 4a) time-dependence R(t), f=2.44 ##########
########### Figure 4b) time-dependence R(t), f=1.7 and R^(-2) plot? ##########


# see magnetopause/run_beta_star_r_mp.sh


########### Figure 5a) time-dependence R(t), f=2.44 ##########
########### Figure 5b) time-dependence R(t), f=1.7 and R^(-2) plot? ##########

# see magnetopause/plot_B_r.sh

########### Figure 6) Summary timeseries of 5 parameters for the north (and south?) pole ##########

# see carrington_plot_timeseries.sh





