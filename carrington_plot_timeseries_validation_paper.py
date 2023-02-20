import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path

import os, sys
from read_vsc_data import *    # read_vsc_data reads txt file output from vlsvintpol.py into a dictionary
import numpy as np

import argparse
import keogram as keogram
from copy import deepcopy
import pandas as pd


# Program purpose: makes timeseries plots of the main inferred parameters indicative of solar storms

parser = argparse.ArgumentParser()
parser.add_argument('-run', default='EGL', help="the Vlasiator run, e.g. egl or EGL" )
parser.add_argument('-n_sc', help="number of coordinates in the data files" )
parser.add_argument('-t', default='t_sec', help="name of time variable" )
args = parser.parse_args()

#import pdb; pdb.set_trace()

nvar = 4

#runs = ['EGL', 'EGI']
runs = ['EGL', 'EGI', 'EGP']
nruns = len(runs)
fig, axes = plt.subplots(nvar, nruns, figsize= (12,7))      #axes is a 1D array with npos elements

print('start')

for c in range(2): # day/night column

    print(c)
    for run in runs:
    #for run in ['EGL', 'EGI']:
    
        root_dir = '/wrk-vakka/users/horakons/carrington/data/'
        save_dir = '{}{}/aurora_plot_data/'.format(root_dir, run.upper())         #CUSTOM
        
        plot_dir = '/wrk-vakka/users/horakons/carrington/plots/{}_validation_paper/'.format(run)
        plot_file = plot_dir  + 'carrington_summary_timeseries_validation_paper.pdf'
        
        mpl.rcParams.update({'font.size': 16}) 
        plt.rcParams['axes.labelsize'] = 16
        plt.rcParams['axes.titlesize'] = 16

        restore_data = False   # MUCH faster if True, re-calculate data if false
        if restore_data==True:   # don't reprocess data, use pickle file instead 
            print('Restoring pre-saved file...')
            tmp = restore(save_dir+'carrington_plot_timeseries.pickle')
            t = tmp['t']
            dct_plot = tmp['dct_plot']
        else:
            datafile = save_dir + 'data.csv'        
            
            if args.n_sc is None:
                n_sc_input = None
            else:
                n_sc_input= int(args.n_sc)

            print('n_sc_input', n_sc_input)
            print('tvar', args.t)
            print('reading')
            dct, nvars = read_vsc_data(datafile, n_sc = n_sc_input, tvar = args.t)
            print('done reading')    

            t = dct[args.t][:,0]
            ntimes = t.size
            
            # TODO: implement this separately for North and South poles
            
            proton_DNF_10keV_max_day = np.zeros([ntimes])
            proton_DNF_10keV_max_night = np.zeros([ntimes])
            #proton_DNF_15811_eV_max = np.zeros([ntimes])
            last_closed_north_deg = np.zeros([ntimes])
            #first_open_north_deg = np.zeros([ntimes])
            first_open_north_dayside_deg = np.zeros([ntimes])
            first_open_north_nightside_deg = np.zeros([ntimes])
            jpar_north_max_day = np.zeros([ntimes])
            jpar_north_max_night = np.zeros([ntimes])
            #last_closed_south_deg = np.zeros([ntimes])
            #first_open_south_deg = np.zeros([ntimes])
            #jpar_south_max = np.zeros([ntimes])
            magnetopause_standoff_RE = np.zeros([ntimes])
            
            # calculate quantities of interest
            
            theta_GSE_deg = dct['theta_GSE_deg'][0,:]
            phi_GSE_deg = dct['phi_GSE_deg'][0,:]
            
            
            #theta_GSE_deg goes from 0 (north pole)--> 180 (south pole)
            #phi_GSE_deg goes from -180 to 180  (0 = noon)
            
            
            # Key names:
            #r_GSE_m,theta_GSE_deg,phi_GSE_deg,r_GSE_traced_m,theta_GSE_traced_deg,phi_GSE_traced_deg,t_sec,open_vs_closed,proton_DNF_00500_eV,proton_DNF_00889_eV,proton_DNF_01581_eV,proton_DNF_02811_eV,proton_DNF_05000_eV,proton_DNF_10keV,proton_DNF_15811_eV,proton_DNF_28117_eV,proton_DNF_50000_eV,$(B_0 / B)J_\parallel$  $[A/km^2]$,Magnetopause_radius_RE
            
            if run == 'EGP':
                keys = ['proton_DNF_09999_eV', 'open_vs_closed', '$(B_0 / B)J_\parallel$  $[A/km^2]$']   #, 'Magnetopause_radius_RE\n']
            else: # EGI, EGL
                keys = ['proton_DNF_10keV', 'open_vs_closed', '$(B_0 / B)J_\parallel$  $[A/km^2]$']   #, 'Magnetopause_radius_RE\n']
            #keys = ['proton_DNF_15811_eV', 'open_vs_closed', '$(B_0 / B)J_\parallel$  $[A/km^2]$']   #, 'Magnetopause_radius_RE\n']
            #keys = ['proton_DNF_15811_eV', 'open_vs_closed', 'Jpar']

            print('phi_GSE min:', np.nanmin(dct['phi_GSE_deg']))
            print('phi_GSE max:', np.nanmax(dct['phi_GSE_deg']))
            

            mask_dayside = (np.abs(phi_GSE_deg) <= 90) & (theta_GSE_deg < 90)  #(north pole)
            mask_nightside = (np.abs(phi_GSE_deg) > 90) & (theta_GSE_deg < 90) 
            #inds_dayside, = np.where( (np.abs(dct['phi_GSE_deg'].flatten()) <= 90) & (dct['theta_GSE_deg'].flatten() < 90) ) #(north pole)
            #inds_nightside, = np.where( (np.abs(dct['phi_GSE_deg'].flatten()) > 90) & (dct['theta_GSE_deg'].flatten() < 90) )
                
            for i in range(ntimes):
                print(i)
                # 1. proton DNF
                #proton_DNF_10keV_max_day[i] = np.nanmax(dct[keys[0]].flatten()[inds_dayside])
                #proton_DNF_10keV_max_night[i] = np.nanmax(dct[keys[0]].flatten()[inds_nightside])
                proton_DNF_10keV_max_day[i] = np.nanmax(dct[keys[0]][i,:][mask_dayside])
                proton_DNF_10keV_max_night[i] = np.nanmax(dct[keys[0]][i,:][mask_nightside])
                #proton_DNF_15811_eV_max[i] = np.nanmax(dct[keys[0]][i,:])
                # 2. last closed field line
                #inds_closed_north, = np.where( (dct[keys[1]][i,:] == 0.0) & (theta_GSE_deg < 90))
                #inds_closed_south, = np.where( (dct[keys[1]][i,:] == 0.0) & (theta_GSE_deg > 90))
                #last_closed_north_deg[i] = 90-np.nanmin(theta_GSE_deg[inds_closed_north])
                #last_closed_south_deg[i] = 90-np.nanmax(theta_GSE_deg[inds_closed_south]) 
                #inds_open_north, = np.where( (dct[keys[1]][i,:] == 1.0) & (theta_GSE_deg < 90))
                inds_open_north_dayside, = np.where( (dct[keys[1]][i,:] == 1.0) & (theta_GSE_deg < 90) & 
                                                     (np.abs(phi_GSE_deg) == np.nanmin(np.abs(phi_GSE_deg))) )
                inds_open_north_nightside, = np.where( (dct[keys[1]][i,:] == 1.0) & (theta_GSE_deg < 90) & 
                                                       ( np.abs(phi_GSE_deg) == np.nanmax(np.abs(phi_GSE_deg)) ) )
                #inds_open_south, = np.where( (dct[keys[1]][i,:] == 1.0) & (theta_GSE_deg > 90))
                #first_open_north_deg[i] = 90-np.nanmax(theta_GSE_deg[inds_open_north])
                first_open_north_dayside_deg[i] = 90-np.nanmax(theta_GSE_deg[inds_open_north_dayside])
                first_open_north_nightside_deg[i] = 90-np.nanmax(theta_GSE_deg[inds_open_north_nightside])     #REIMPLEMENT THIS!!
                #first_open_south_deg[i] = 90-np.nanmin(theta_GSE_deg[inds_open_south]) 
                # 3. parallel current
                #jpar_north_max_day[i] = np.nanmax( np.abs(dct[keys[2]].flatten()[inds_dayside]) )
                #jpar_north_max_night[i] = np.nanmax( np.abs(dct[keys[2]].flatten()[inds_nightside]) )
                jpar_north_max_day[i] = np.nanmax( np.abs(dct[keys[2]][i,:][mask_dayside]) )
                jpar_north_max_night[i] = np.nanmax( np.abs(dct[keys[2]][i,:][mask_nightside]) )
                #jpar_south_max[i] = np.nanmax( np.abs(dct[keys[2]][i,mask_south]) )
                # 4. magnetopause standoff
                #magnetopause_standoff_RE[i] = np.nanmin( np.abs(dct[keys[3]][i,:]) ) 
            dct_plot = {}
            dct_plot['DNF'] = {'variables_day':[proton_DNF_10keV_max_day], 'variables_night':[proton_DNF_10keV_max_night],
                               'labels':['8.9 keV'],   # 8891 eV
            #dct_plot['DNF'] = {'variables_day':[proton_DNF_15811_eV_max], 'labels':['15811 eV'],
                               'ylabel':'max. DNF\n$[cm^{-2}s^{-1}$\n$sr^{-1}eV^{-1}]$', 'ylim':[0, 1500] }  # colors = ['r']
                               #'ylabel':r'max. DNF $[cm^{-2}s^{-1}sr^{-1}eV^{-1}]$' }  # colors = ['r']
            #dct_plot['Last closed lat. (deg.)'] = last_closed_north_deg
            #dct_plot['First open lat. (deg.)'] = first_open_north_deg
            dct_plot['OCB'] = {'variables_day':[first_open_north_dayside_deg], 'variables_night':[first_open_north_nightside_deg],
                               'labels':['noon meridian'], 'ylabel':'OCB\n[deg.]', 'ylim':[68, 78] }
            #dct_plot['OCB'] = {'variables_day':[first_open_north_dayside_deg, first_open_north_nightside_deg], 
                               #'labels':['noon','midnight'], 'ylabel':'OCB [deg.]', 'ylim':[68, 78] }
            dct_plot['FAC'] = {'variables_day':[jpar_north_max_day],'variables_night':[jpar_north_max_night], 'labels':[None], 
                               'ylabel':'max. FAC\n$[A km^{-2}]$', 'ylim':[0.5, 2.5]}  # colors = ['r']
                               #'ylabel':r'$(B_0 / B) J_\parallel [A km^{-2}]$'}  # colors = ['r']
            #dct_plot['Magnetopause standoff (RE)'] = magnetopause_standoff_RE
       

        dct_plot['DNF']['ylabel'] = 'max. DNF\n$[cm^{-2}s^{-1}$\n$sr^{-1}eV^{-1}]$'
        dct_plot['OCB']['ylabel'] = 'OCB\n[deg.]'
        dct_plot['FAC']['ylabel'] = 'max. FAC\n$[A km^{-2}]$'
        
        t0 = 857   # nominal pulse arrival time, EGL
        
        # vv         add Emilia's plot here        vv
        df_sc = pd.read_csv('/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/mp_nose_x_new_'+run+'_backup.csv')   # data generated with updated version of beta_star_r_mp_new.py
        #df_sc = pd.read_csv('/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/mp_nose_x.csv')   # data generated with updated version of beta_star_r_mp.py
        dct_plot['Mpos'] = {'variables_day':[np.array(df_sc['x'])], 'variables_night':[np.array(df_sc['x'])*0], 'labels':[None], 'ylabel':'R\n[$R_E$]',  'ylim':[8, 11]}

        print(dct_plot['Mpos']['variables_day'][0])
        print(dct_plot['Mpos']['variables_night'][0])
        print(dct_plot['Mpos']['labels'][0])
        print(dct_plot['Mpos']['labels'][0])
        print('saving')
        #SAVE data
        if restore_data == False:
            save(save_dir+'carrington_plot_timeseries.pickle', t=t, dct_plot=dct_plot)
        
        #dct_plot_keys = dct_plot.keys()
        dct_plot_keys = ['OCB', 'FAC', 'DNF', 'Mpos']
        if c==0:
            figletter = ['(a)', '(b)', '(c)', '(d)']
        else:
            figletter = ['(e)', '(f)', '(g)', '']

        
        #titles = ['DNF, 15811 eV', 'Field line topology', '$(B_0 / B)J_\parallel$  $[A/km^2]$']
        for i, key in enumerate(dct_plot_keys):
            print(key)
            try:
                if figletter[i] != '':
                    if key == 'OCB':
                        if c==0:
                            axes[i,c].set_title('Dayside', fontsize = 18)
                        else:
                            axes[i,c].set_title('Nightside', fontsize = 18)
                    if c==0:
                        for j in range(len( dct_plot[key]['variables_day'] )):
                            print('test:')
                            print(type(dct_plot[key]['variables_day'][j]))
                            print(type(dct_plot[key]['labels'][j]))
                            print(j)
                            print('shapes:', t.shape, dct_plot[key]['variables_day'][j].shape)
                            y = dct_plot[key]['variables_day'][j]
                            print('SIZES:', t.size, ' ', y.size)
                            if y.size < t.size:
                                axes[i,c].plot(t[0:y.size-1], dct_plot[key]['variables_day'][j], label = dct_plot[key]['labels'][j])
                            else:
                                axes[i,c].plot(t, dct_plot[key]['variables_day'][j], label = dct_plot[key]['labels'][j])
                            print('DONE plotting day')
                    else:
                        for j in range(len( dct_plot[key]['variables_night'] )):
                            print(j)
                            axes[i,c].plot(t, dct_plot[key]['variables_night'][j], label = dct_plot[key]['labels'][j])
                            print('DONE plotting night')
                    #if dct_plot[key]['labels'][0] is not None:
                    #    axes[i,c].legend(framealpha=0.5, loc='center left')
                    axes[i,c].set_xlim([600,1800])
                    axes[i,c].set_xbound([600,1800])
                    axes[i,c].set_ylim(dct_plot[key]['ylim'])
                    axes[i,c].plot([t0, t0], dct_plot[key]['ylim'], linestyle = ':', color = 'green')
                    if c == 0: # left column
                        if i <= nvar-2:
                            #axes[i,c].set_visible(False)
                            axes[i,c].set_xlabel('')
                            for xlabel_i in axes[i,c].get_xticklabels():
                                xlabel_i.set_visible(False)
                        else:
                            axes[i,c].set_xlabel('time [seconds]', fontsize = 16)
                        axes[i,c].set_ylabel(dct_plot[key]['ylabel'], fontsize = 16)
                    else:  # right column
                        if i <= nvar-3:
                            #axes[i,c].set_visible(False)
                            axes[i,c].set_xlabel('')
                            for xlabel_i in axes[i,c].get_xticklabels():
                                xlabel_i.set_visible(False)
                        else:
                            axes[i,c].set_xlabel('time [seconds]', fontsize = 16)
                        axes[i,c].set_ylabel('')
                        for ylabel_i in axes[i,c].get_yticklabels():
                            ylabel_i.set_visible(False)
                    extraticks = [500, 621, t0, 1000, 1500, 1760]
                    axes[i,c].locator_params(tight=True, nbins=4)
                    axes[i,c].set_xticks(list(axes[i,c].get_xticks()) + extraticks)
                    axes[i,c].grid(axis='y', linewidth=1)  # , linestyle = 
                    ax2 = axes[i,c].twinx()
                    color = 'black'
                    ax2.set_ylabel(figletter[i], color = color)
                    plt.tick_params(right = False , labelright = False)
                    plt.subplots_adjust(hspace=0.8, wspace=0.6, left = 0.25)
                    #axes[i,c].annotate(figletter[i], 0.1, 0.9, xycoords = 'figure fraction') doesn't plot?
                    #fig.tight_layout()
                else:
                    axes[i,c].axis('off')
                    from matplotlib.patches import Patch
                    from matplotlib.lines import Line2D
                    legend_elements = [Line2D([0], [0], color='C0', lw=4, label='Pulse'),
                                       Line2D([0], [0], color='C1', lw=4, label='Control'), 
                                       Line2D([0], [0], color='C2', lw=4, label='EGP')]
                    axes[i,c].legend(handles=legend_elements, loc='center')
            except:
                print('key {} error!'.format(key))
plt.tight_layout()
axes[0,0].text(800, 6.7, r'$t_0$')
axes[0,1].text(800, 25.5, r'$t_0$')
print(plot_file)
mkdir_path(plot_file)
plt.savefig(plot_file)
plt.close()













