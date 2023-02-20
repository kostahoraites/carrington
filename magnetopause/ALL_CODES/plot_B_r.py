import pytools as pt
import numpy as np
import pandas as pd
from myutils import get_vlsvfile_fullpath
import matplotlib.pyplot as plt
import csv
from scipy.signal import savgol_filter
from tsyganenko import tsyganenko_b
from xibi import get_b_xaxis


# v IMPLEMENT THIS v
#
#before and after parameters
fileIndexs = [621, 1760]    #857
Dsts = [-33, -33]
Bz_imfs = [-5., -10.]
N_sws = [1., 4.]
rmags = [10.36, 8.22]    #8.221087695799728,    10.357921440969816
run = 'EGL'


#for fileIndex, Dst, Bz_imf, N_sw in zip(fileIndexs, Dsts, Bz_imfs, N_sws):
# [...]
# indent and loop
# input  Dst, Bz_imf, N_sw  into call to tyganenko_b
# construct files xibi_{fileIndex}.csv for different times
# B_r = open('xibi.csv') # change this line (below)
# change the filenames in savefig() calls 
# add "t={fileIndex} sec" to titles / legends
# legends show tsyganenko and Vlasiator
# 2 separate plots (pdf/eps), put them side by side in draft
#
# ^               ^


for fileIndex, Dst, Bz_imf, N_sw, rmag in zip(fileIndexs, Dsts, Bz_imfs, N_sws, rmags):

    B_dict = get_b_xaxis(run, fileIndex, save = False)
    x = B_dict['x']
    B = B_dict['B']
    
    #READ DATA LATER:
    #B_r = open('xibi.csv')
    #B_r = open('xibi_ave.csv')
    #csvreader = csv.reader(B_r)
    #header = []
    #header = next(csvreader)
    
    #x = []
    #B = []
    #for row in csvreader:
    #    x.append(float(row[0]))
    #    B.append(float(row[1]))
    
    B_log10 = np.log10(np.array(B))
    
    #calculate Tsyganenko B
    xarr = np.array(x)
    yarr = xarr*0
    zarr = xarr*0
    bx_tsyg, by_tsyg, bz_tsyg = tsyganenko_b(xarr, yarr, zarr, Txx = 't01', InternalB='dipole', Dst = Dst, Kp = 4, Vx_sw = -750., N_sw = N_sw, Bx_imf = 0., By_imf = 0., Bz_imf = Bz_imf )
    
    B_tsyg = np.sqrt(bx_tsyg**2 + by_tsyg**2 + bz_tsyg**2) * 1e-9    # [T]
    
    lnr_tsyg = np.log(xarr)
    lnB_tsyg = np.log(B_tsyg)
    B_log10_tsyg = np.log10(np.array(B_tsyg))
    
    nx = xarr.size
    dlnB_dlnr_tsyg = (lnB_tsyg[1:nx] - lnB_tsyg[0:nx-1]) / (lnr_tsyg[1:nx] - lnr_tsyg[0:nx-1]) 
    
    print(dlnB_dlnr_tsyg)
    print(xarr.size)
    print(dlnB_dlnr_tsyg.size)
    
    #SMOOTH B:
    #B = savgol_filter(np.array(B), 9, 3) #these numbers gave the best result!
    
    i = 0
    dB_list = []
    dx_list = []
    while i+1 < len(B):
        dB = B[i+1] - B[i]
        dB_list.append(dB)
        dx = x[i+1] - x[i]
        dx_list.append(dx)
        i += 1
    
    i = 0
    B_new = []
    x_new = []
    while i+1 < len(x):
        x_n = (x[i+1] + x[i])/2
        x_new.append(x_n)
        B_n = (B[i+1] + B[i])/2
        B_new.append(B_n)
        i += 1
    
    i = 0
    slope_list = []
    while i < len(dB_list):
        slope = x_new[i]/B_new[i]*dB_list[i]/dx_list[i]
        slope_list.append(slope)
        i += 1

    plt.rcParams.update({'font.size': 14})   
    plt.figure(figsize=(6,7))
    #plt.subplot(2,1,1)
    #plt.plot([5.4, 12.2], [-2.4, 0], 'k--', linewidth=1.5)
    plt.plot(x_new, dlnB_dlnr_tsyg, label = 'Tsyganenko (2001)')
    plt.scatter(x_new, slope_list, color = 'r', marker='.', label = 'Vlasiator')
    plt.xscale('log')
    plt.xlabel(r'$r\: [R_E]$', fontsize='18')
    plt.ylabel(r'index $(d \ln B / d \ln r)$', fontsize='18')
    plt.title('B scaling index, t={} sec'.format(fileIndex))
    plt.axvline(x=rmag, color='orange', ls='--', label=r'Magnetopause')
    #plt.axvline(x=8.169, color='r', ls='--', label=r'Magnetopause')
    plt.axvline(x=4.7, color='y', ls='--', label=r'$r_{min}$')
    plt.grid(axis='y', linewidth=1.5, linestyle=':')
    plt.xlim(1e0, 5e1)
    plt.ylim(-4, 2)
    plt.legend()
    plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/B_r_slope_{}.pdf'.format(fileIndex))
    plt.close()

    plt.figure(figsize=(6,7))
    #plt.subplot(2,1,2)
    plt.plot(x, B_log10_tsyg, label = 'Tsyganenko (2001)')
    plt.scatter(x, B_log10, color = 'r', marker='.', label = 'Vlasiator')
    plt.axvline(x=rmag, color='orange', ls='--', label=r'Magnetopause')
    #plt.axvline(x=8.169, color='r', ls='--', label=r'Magnetopause')
    plt.axvline(x=4.7, color='y', ls='--', label=r'$r_{min}$')
    plt.plot([1, 10], [-3, -6], color = 'black', linewidth=1.5)
    plt.text(3, -3.9, r'$r^{-3}$')
    plt.title('B(r), t={} sec'.format(fileIndex))
    plt.xscale('log')
    plt.xlabel(r'$r\: [R_E]$', fontsize='18')
    plt.ylabel(r'$\log_{10} (B\: [T])$', fontsize='18')
    plt.xlim(1e0, 5e1)
    plt.legend()
   
    plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/B_r_{}.pdf'.format(fileIndex))
    plt.close()

