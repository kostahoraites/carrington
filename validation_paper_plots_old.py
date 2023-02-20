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
import pandas
from time import time
from carrington import *
from tsyganenko import *


mpl.rcParams.update({'font.size': 18})   # increase font size

#plot_dir = '/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/'


# example call from bash (see carrington.sh):
# python carrington.py -run EGI -var 1 4 6

# Input parameters
parser = argparse.ArgumentParser()

parser.add_argument('-nproc', default=1, help="number of processors to use " )
global ARGS
ARGS = parser.parse_args()
#nproc = ARGS.nproc



########### Figure 1a) 3D perspective plot of density ##########



#EXAMPLE 5: plot_threeslice

import pytools as pt

# Defining source and output file locations
run = 'EGL'
dim = '3D'
bulk = 'bulk1.egl'
bulkLocation = '/wrk/group/spacephysics/vlasiator/{}/{}/bulk/'.format(dim, run)
fluxLocation = '/wrk/group/spacephysics/vlasiator/{}/{}/flux/'.format(dim, run)
outputLocation = '/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/'
j=857                       #EGL --- pressure pulse arrival time
bulkname = '{}.{}.vlsv'.format(bulk,str(j).zfill(7))

#pt.plot.plot_threeslice(filename=bulkLocation+bulkname, var ='proton/vg_rho', colormap = 'plasma', vmin = 1e4,vmax=5e6, step=j, outputdir=outputLocation,outputfile='test_threeslice_{}_{}.jpg'.format(run,j), Earth=1,cutpointre=[0,0,0], slices='yz' )
pt.plot.plot_threeslice(filename=bulkLocation+bulkname, var ='vg_pdyn', colormap = 'plasma', vmin = 1e-10,vmax=5e-9, step=j, outputdir=outputLocation,outputfile='test_threeslice_pdyn_{}_{}.jpg'.format(run,j), Earth=1,cutpointre=[0,0,0], slices='yz' )

plt.close()



########## Figure 1b) Pearson correlation, best fit with OMNI ###########


########## Figure 2 (setup) ###########

run = 'EGL'
global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations
global plot_dir
plot_dir = '/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/'
global mu_0
mu_0 = 4e-7 * np.pi




def f_ocb_before(phi_i):
    return tsyganenko_ocb(phi_i, lat_range = [60,80], nsteps = 10,
                          Txx = 't01', Dst = -30, Kp = 4, Vx_sw = -750., N_sw = 1., Bx_imf = 0., By_imf = 0., Bz_imf = -5.,
                          R_inner = 1., R_outer=15., dir = None, maxloop = 10000)


def f_ocb_after(phi_i):
    return tsyganenko_ocb(phi_i, lat_range = [60,80], nsteps = 10,
                          Txx = 't01', Dst = -80, Kp = 4, Vx_sw = -750., N_sw = 4., Bx_imf = 0., By_imf = 0., Bz_imf = -10.,
                          R_inner = 1., R_outer=15., dir = None, maxloop = 10000)


# define quantities normally contained in ARGS
    
phimin = -180.
phimax = 180.
latmin = 60.  # (CHANGE THIS!)
latmax = 90.
nlat = 121   #121
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
#ocb_before = pool.map(f_ocb_before, phi_1d)
#pool.close()
#pool.join()

#pool = Pool(int(ARGS.nproc))
#ocb_after = pool.map(f_ocb_after, phi_1d)
#pool.close()
#pool.join()

#np.savetxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_before.txt', ocb_before)    
#np.savetxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_after.txt', ocb_after)

## ^        ^

ocb_tsyg_before = np.loadtxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_before.txt')    
ocb_tsyg_after = np.loadtxt('/wrk-vakka/users/horakons/carrington/data/sidecar/ocb_EGL_after.txt')
phi_1d_tsyg = np.linspace(-np.pi, np.pi, num = ocb_tsyg_after.size)
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
        for j in range(1, nlat):
            if latcut[j] != latcut[j-1]:  # find topology changes, closed <--> open 
                 ocb_1d[i] = ( (lat_1d[j] + lat_1d[j-1]) / 2.)  # already in degrees
    return ocb_1d

data = {'before':{'fileIndex':621, 'ocb_tsyg':ocb_tsyg_before}, 'after':{'fileIndex':1760, 'ocb_tsyg':ocb_tsyg_after}}


labels = ['a)', 'c)', 'e)', 'b)', 'd)', 'f)']
ilabel = 0


for key in data: # 'before', 'after'
    
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
    
    #trace_coord = sidecar(trace_coord)        # SIDECAR: read function output from hard disk if calulated before
    
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
    ax.set_title('OCB, t={} sec'.format(fileIndex))
    plt.text(0.1, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 20)
    ilabel += 1
    #plt.legend()
    plt.savefig(plot_dir+'ocb_{}_azim.pdf'.format(key))
    plt.close()

    data[key]['ocb_2d_tf'] = fOO.data_list[0]     # boolean, [nlat, nphi]

    data[key]['ocb_1d'] = ocb_2d_to_1d(data[key]['ocb_2d_tf'], lat_1d)
    data[key]['ocb_1d'] = ocb_2d_to_1d(data[key]['ocb_2d_tf'], lat_1d)

    ########### Figure 2b) FACs ##########

    # Field aligned currents (FACS)

    JO_traced = FAC_Obj(run, f, coord_list_aurora, fg_b = fg_b, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data = False)
    ax = JO_traced.make_plot_data( plot_data=True, plots = [1], pltclose=False)  # polar plot
    #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
    ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label=ocb_label, color = "green" )
    plt.text(0.1, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 20)
    ilabel += 1
    ax.set_title('FAC, t={} sec'.format(fileIndex))
    #plt.legend()
    plt.savefig(plot_dir+'fac_{}_azim.pdf'.format(key))
    plt.close()

    ########### Figure 2c) precipitation ##########

    # proton differential energy flux

    pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data=False)
    i = 5
    ax = pFO_traced.make_plot_data( plot_data=True, plots = [i], pltclose=False)  # polar plot
    #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
    ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label = ocb_label, color = "green" )
    plt.text(0.1, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 20)
    ilabel += 1
    proton_energy = f.read_parameter('proton_PrecipitationCentreEnergy{}'.format(i))
    ax.set_title('DNF, {:.1f} keV, t={} sec'.format(proton_energy/1000., fileIndex))
    #plt.legend()  '{:.1f}'.format(8891 / 1000.)
    plt.savefig(plot_dir+'def_{}_{}_azim.pdf'.format(key, int(proton_energy)))
    plt.close()
    
    pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = False, plot_data=False)
    i = 6
    ax = pFO_traced.make_plot_data( plot_data=True, plots = [i], pltclose=False)  # polar plot
    #ax.plot(phi_azim_tsyg_plot, 90-data[key]['ocb_tsyg'], color="black", linewidth=2, label = tsyg_label)
    ax.plot(phi_azim_plot, 90-data[key]['ocb_1d'], label = ocb_label, color = "green" )
    plt.text(0.1, 0.9, labels[ilabel], transform=plt.gcf().transFigure, fontsize = 20)
    proton_energy = f.read_parameter('proton_PrecipitationCentreEnergy{}'.format(i))
    ax.set_title('DNF, {:.1f} keV, t={} sec'.format(proton_energy/1000., fileIndex))
    #plt.legend()
    plt.savefig(plot_dir+'def_{}_{}_azim.pdf'.format(key, int(proton_energy)))
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


#RESULT: data has 'fileIndex', 'ocb_2d_tf', and 'fg_b' keys



ocb_1d_before = data['before']['ocb_1d']
ocb_1d_after = data['after']['ocb_1d']

#plt.plot( phi_azim_plot, ocb_1d_before)
#plt.savefig(plot_dir+'ocb_before_line.png')
#plt.close()

#plt.plot( phi_azim_plot, ocb_1d_after)
#plt.savefig(plot_dir+'ocb_after_line.png')
#plt.close()

fig, ax = plt.subplots()
plt.xlabel(r'MLT [hours]')
plt.ylabel(r'OCB Latitude [$^{\circ}$]')
plt.title(r'OCB vs. MLT')
plt.plot(phi_mlt, ocb_1d_before, label='OCB (before)' )
plt.plot(phi_mlt, ocb_1d_after, label='OCB (after)' )
ax.set_xticks([6,8,10,12,14,16,18])
ax.set_xlim([6,18])
#ax.set_xticks([0,4,8,12,16,20,24])
#plt.plot(phi_mlt_tsyg, ocb_tsyg_before, label='Tsyg. (before)' )
#plt.plot(phi_mlt_tsyg, ocb_tsyg_after, label='Tsyg. (after)' )
#plt.xlim([0,24])
filename = plot_dir + 'ocb_vs_mlt.pdf'
plt.legend()
plt.tight_layout() 
plt.savefig(filename)
plt.close()





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





