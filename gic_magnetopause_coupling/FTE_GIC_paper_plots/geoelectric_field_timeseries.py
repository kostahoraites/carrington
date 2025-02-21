'''
Before running cell blocks below, requires running biot_savart.py 
(option #1 in main()) to generate total ground magnetic field files
'''

import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import ftest as ft
import pytools as pt
from myutils import cartesian_to_spherical_vector, cartesian_to_spherical, spherical_to_cartesian, mkdir_path, timer, get_vlsvfile_fullpath

R_EARTH = 6371000.

run = "FHA"  # FHA, FIA, EGL
if run == "FIA":
    dir = "/wrk-vakka/group/spacephysics/vlasiator/3D/FIA/bulk_sidecars/ig_B"
elif run == "FHA":
    dir = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/ig_B"
elif run == "EGL":
    dir = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/sidecars/ig_B"
f = ft.f("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/ig_B/ionosphere_B_sidecar_FHA.0000784.vlsv")

def E_horizontal(dB_dt, pos, time, sigma = 1e-3, method = 'liu'):
    '''
        Calculate the horizontal electric field by integrating components of dB/dt
            References: Cagniard et al 1952 (eq. 12), Pulkinnen et al 2006 (eq. 19)
        Inputs:
            dB_dt: cartesian dB/dt    [T/s] array dimension [3, len(time)]
            pos: cartesian position [m] 1D 3-element array, vector position
            time: 1D array of times [s], monotonically increasing

        Keywords:
            sigma = ground conductivity (siemens/meter)
            method:
                'liu': use integration method described in Liu et al., (2009) doi:10.1029/2008SW000439, 2009
                       this method is exact for piecewise linear B (i.e., piecewise constant dB/dt)
                'RH-riemann': use right-handed Riemann sum.
    '''
    mu_0 = 1.25663706e-6    # permeability of free space
    E_north = np.zeros(time.size)
    E_east = np.zeros(time.size)
    dB_dt_r, dB_dt_theta, dB_dt_phi = cartesian_to_spherical_vector(dB_dt[0,:], dB_dt[1,:], dB_dt[2,:], pos[0], pos[1], pos[2])
    dB_dt_north = -dB_dt_theta; dB_dt_east = dB_dt_phi
    t0 = time[1] - time[0]
    for i in range(0, time.size):
        t = time[i]   # monotonically increasing from t[0]
        tp = time[1:i+1]
        dt = tp - time[0:i]
        if method == 'liu':  # implement Liu et al (2009), eq. 5
            # Possible upgrade to the integration would be to fit the dB/dt time-series to a higher-order polynomial
            # and integrate that exactly over the finite time domain
            t_1 = t - tp
            t_2 = t - time[0:i]   # elementwise, t_2 > t_1
            dB_north = dB_dt_north[0:i] * dt[0:i]
            dB_east = dB_dt_east[0:i] * dt[0:i]
            E_north[i] = np.sum(-(2. / np.sqrt(np.pi * mu_0 * sigma * dt[0:i])) * dB_east * (np.sqrt(t_2) - np.sqrt(t_1) ) ) 
            E_east[i] = np.sum((2. / np.sqrt(np.pi * mu_0 * sigma * dt[0:i])) * dB_north * (np.sqrt(t_2) - np.sqrt(t_1) ) ) 
        elif method == 'RH-riemann':
            if i != 0:
                E_north[i] = -(1. / np.sqrt(np.pi * mu_0 * sigma)) * np.sum(dt * dB_dt_east[1:i+1] / np.sqrt(t-tp + t0))
                E_east[i] = (1. / np.sqrt(np.pi * mu_0 * sigma)) * np.sum(dt * dB_dt_north[1:i+1] / np.sqrt(t-tp + t0))  # note the sign
    return E_north, E_east



if run == "FHA":
    nmin = 501    # 1001
    nmax = 1612    # 1100
elif run == "FIA":
    nmin = 1
    nmax = 817         #865 (files 818-819 missing)
elif run == "EGL":
    nmin = 621
    nmax = 1760
time = np.linspace(nmin, nmax, nmax - nmin + 1)

pos = f.read_variable('ig_r')   # ionospheric grid.  array dimensions (43132, 3)
#ig_B_arr = np.ndarray([pos.shape[0], pos.shape[1], nmax - nmin + 1])
ig_dB_dt_arr = np.ndarray([pos.shape[0], pos.shape[1], nmax - nmin + 1])
ig_B_ionosphere_arr = ig_dB_dt_arr * 0.
ig_dB_dt_ionosphere_arr = ig_dB_dt_arr * 0.
ig_B_inner_arr = ig_dB_dt_arr * 0.
ig_dB_dt_inner_arr = ig_dB_dt_arr * 0.
ig_B_outer_arr = ig_dB_dt_arr * 0.
ig_dB_dt_outer_arr = ig_dB_dt_arr * 0.

E_north_arr = np.ndarray([pos.shape[0], nmax - nmin + 1])
E_east_arr = np.ndarray([pos.shape[0], nmax - nmin + 1])

#populate B arrays
for i in range(nmin, nmax+1):
    #print(i)
    f = ft.f(dir + "/ionosphere_B_sidecar_{}.{}.vlsv".format(run, str(i).zfill(7)))     # FHA: file indices 501 - 1612
    try:
        ig_B_ionosphere = f.read_variable('ig_B_ionosphere')
        ig_B_ionosphere_arr[:,:,i-nmin] = ig_B_ionosphere
    except:
        print("couldn't read ionospheric data") # for runs without an ionosphere, leave as zeros
    ig_B_inner = f.read_variable('ig_B_inner')
    ig_B_inner_arr[:,:,i-nmin] = ig_B_inner
    ig_B_outer = f.read_variable('ig_B_outer')
    ig_B_outer_arr[:,:,i-nmin] = ig_B_outer

# interpolate across zeros in B arrays (missing data points, something funny about FHA and FIA makes ionosphere write operation unreliable?)
#for arr in [ig_B_arr, ig_B_ionosphere_arr, ig_B_inner_arr, ig_B_outer_arr]:

for arr in [ig_B_ionosphere_arr]:
    try:
        #reduced_arr = np.sum(arr, 0)
        #ind = np.where(reduced_arr[0,:] != 0)[0]
        ind = np.where(arr[0,0,:] != 0)[0]
        print("{} points removed".format(arr.shape[2] - ind.size))
        interp_arr = arr[:,:, ind]  # only keep the non-zero times to conduct the interpolation
        for i in range(arr.shape[0]): # positions
            for j in range(3): # vector components
                arr[i, j, :] = np.interp(time, time[ind], interp_arr[i, j, :], left=None, right=None, period=None)
    except:
        print("error with interpolation. zeroing out array...")

ig_B_arr =  ig_B_ionosphere_arr + ig_B_inner_arr + ig_B_outer_arr

# calculate dB/dt and populate corresponding arrays
for i in range(nmin+1, nmax+1):
    #next compute derivatives
    ig_dB_dt_arr[:,:,i-nmin] = ig_B_arr[:,:,i-nmin] - ig_B_arr[:,:,i-nmin-1]
    ig_dB_dt_ionosphere_arr[:,:,i-nmin] = ig_B_ionosphere_arr[:,:,i-nmin] - ig_B_ionosphere_arr[:,:,i-nmin-1]
    ig_dB_dt_inner_arr[:,:,i-nmin] = ig_B_inner_arr[:,:,i-nmin] - ig_B_inner_arr[:,:,i-nmin-1]
    ig_dB_dt_outer_arr[:,:,i-nmin] = ig_B_outer_arr[:,:,i-nmin] - ig_B_outer_arr[:,:,i-nmin-1]

'''
#i_pos = 0
#E_north, E_east = E_horizontal(ig_dB_dt_arr[i_pos,:,:], pos[i_pos,:], time, sigma = 1e-3, method = 'liu')
for i_pos in range(ig_dB_dt_arr.shape[0]):
    E_north, E_east = E_horizontal(ig_dB_dt_arr[i_pos,:,:], pos[i_pos,:], time, sigma = 1e-3, method = 'liu')
    E_north_arr[i_pos,:] = E_north
    E_east_arr[i_pos,:] = E_east
    print(i_pos)

# write geoelectric field to .vlsv

save_dir = '/wrk-vakka/users/horakons/carrington/data/{}/GIC/'.format(run)
f_iono = pt.vlsvfile.VlsvReader( '/wrk-vakka/group/spacephysics/vlasiator/temp/ionogrid_FHA.vlsv' )

for i, t in enumerate(time):
    # write to file
    filename_vlsv = save_dir + 'ionosphere_gic_{}_{}.vlsv'.format(run, str(int(t)).zfill(7))
    mkdir_path(filename_vlsv)
    writer = pt.vlsvfile.VlsvWriter(f_iono, filename_vlsv)
    writer.write(pos,'ig_r','VARIABLE','ionosphere')
    writer.write(E_north_arr[:,i],'ig_E_north','VARIABLE','ionosphere')
    writer.write(E_east_arr[:,i],'ig_E_east','VARIABLE','ionosphere')
'''


# Plot timeseries of the geoelectric field at different latitudes

plt.rcParams["figure.figsize"] = (10, 6)


f = ft.f(get_vlsvfile_fullpath('FHA', 1165))

'''
# calculate geoelectric field at footpoints of Figure 1 flux ropes
#magnetospheric coord1, coord2
coord1 = np.array([10.450000, -1.430000, -2.600000]) * R_EARTH    # 'FLUX ROPE 1'
coord2 = np.array([10.500000, 1.470000, -2.560000]) * R_EARTH     # 'FLUX ROPE 2'
#calculate ionospheric downmapped coordinates
xyz1 = f.read_interpolated_variable('vg_connection_coordinates_bw', coord1)/R_EARTH
xyz2 = f.read_interpolated_variable('vg_connection_coordinates_bw', coord2)/R_EARTH
r1, theta1, phi1 = cartesian_to_spherical( *xyz1)
r2, theta2, phi2 = cartesian_to_spherical( *xyz2)
theta = np.array([theta1, theta2])
phi = np.array([phi1, phi2])
'''
lat_deg = np.array([-79., -80., -81., 79., 80., 81.])
lat = lat_deg * np.pi / 180.   # (auroral) latitude
theta = (np.pi / 2) - lat  # co-latitude
phi = lat * 0 + 0.00000000000001              # noon, klug to make phi positive


theta_deg = theta * 180. / np.pi 
lat_deg = 90. - theta_deg
phi_deg = phi * 180. / np.pi

x0, y0, z0 = spherical_to_cartesian(R_EARTH, theta, phi)

print(x0.size)
for i in range(x0.size):

    # Find nearest neighbor of the ionosphere grid, index by 'ind_min', to the specified lat and phi
    dist = np.sqrt((x0[i] - pos[:,0])**2 + (y0[i] - pos[:,1])**2 + (z0[i] - pos[:,2])**2) 

    ind_min = np.argmin(dist)
    dB_dt_r, dB_dt_theta, dB_dt_phi = cartesian_to_spherical_vector(ig_dB_dt_arr[ind_min, 0,:], ig_dB_dt_arr[ind_min,1,:], ig_dB_dt_arr[ind_min, 2,:], pos[ind_min, 0], pos[ind_min,1], pos[ind_min, 2])
    dB_dt_north = -dB_dt_theta; dB_dt_east = dB_dt_phi

    E_north, E_east = E_horizontal(ig_dB_dt_arr[ind_min,:,:], pos[ind_min,:], time, sigma = 1e-3, method = 'liu')

    '''
    # vv   OLD PLOTS   vv
    print(z0[i])
    print(pos[ind_min,2])
    print(lat_deg[i], np.nanmax(1e6 * E_north_arr[ind_min,:]))
    # PLOT:
    # geoelectic field
    plt.title('GIC Lat. = {} deg., noon, {}'.format(int(lat_deg[i]), run))
    plt.xlabel('time [sec]')
    plt.ylabel(r'Geoelectric field [$\mu$V/m]')
    plt.plot(time, 1e6 * E_north_arr[ind_min,:], label = r'northward E [$\mu$V/m]')
    plt.plot(time, 1e6 * E_east_arr[ind_min,:], label = r'eastward E [$\mu$V/m]')
    plt.ylim([-400, 400])
    plt.legend()
    filename =  '/wrk-vakka/users/horakons/carrington/plots/{}/GIC/geolectric_E_timeseries_lat_{}_{}'.format(run,int(lat_deg[i]),run)
    mkdir_path(filename)
    plt.savefig(filename)
    plt.close()
    # dB/dt
    plt.title('dB/dt Lat. = {} deg., noon, {}'.format(int(lat_deg[i]), run))
    plt.xlabel('time [sec]')
    plt.ylabel(r'Ground magnetic field [nT/s]]')
    plt.plot(time, 1e9 * dB_dt_north, label = r'northward dB/dt [nT/s]')
    plt.plot(time, 1e9 * dB_dt_east, label = r'eastward dB/dt [nT/s]')
    plt.ylim([-8, 8])
    plt.legend()
    filename =  '/wrk-vakka/users/horakons/carrington/plots/{}/GIC/dB_dt_timeseries_lat_{}_{}'.format(run,int(lat_deg[i]),run)
    mkdir_path(filename)
    plt.savefig(filename)
    plt.close()
    # |dB/dt| (components)
    try: # won't work for EGL?
        plt.title('dB/dt Lat. = {} deg., noon, {}'.format(int(lat_deg[i]), run))
        plt.xlabel('time [sec]')
        plt.ylabel(r'Ground magnetic field [nT/s]]')
        dB_dt_r_ionosphere, dB_dt_theta_ionosphere, dB_dt_phi_ionosphere = cartesian_to_spherical_vector(ig_dB_dt_ionosphere_arr[ind_min, 0,:], ig_dB_dt_arr[ind_min,1,:], ig_dB_dt_arr[ind_min, 2,:], pos[ind_min, 0], pos[ind_min,1], pos[ind_min, 2])
        dB_dt_north_ionosphere = -dB_dt_theta_ionosphere; dB_dt_east_ionosphere = dB_dt_phi_ionosphere
        dB_dt_r_inner, dB_dt_theta_inner, dB_dt_phi_inner = cartesian_to_spherical_vector(ig_dB_dt_inner_arr[ind_min, 0,:], ig_dB_dt_inner_arr[ind_min,1,:], ig_dB_dt_inner_arr[ind_min, 2,:], pos[ind_min, 0], pos[ind_min,1], pos[ind_min, 2])
        dB_dt_north_inner = -dB_dt_theta_inner; dB_dt_east_inner = dB_dt_phi_inner
        dB_dt_r_outer, dB_dt_theta_outer, dB_dt_phi_outer = cartesian_to_spherical_vector(ig_dB_dt_outer_arr[ind_min, 0,:], ig_dB_dt_outer_arr[ind_min,1,:], ig_dB_dt_outer_arr[ind_min, 2,:], pos[ind_min, 0], pos[ind_min,1], pos[ind_min, 2])
        dB_dt_north_outer = -dB_dt_theta_outer; dB_dt_east_outer = dB_dt_phi_outer
        plt.plot(time, np.abs(1e9 * np.sqrt(dB_dt_north**2 + dB_dt_east**2) ), label = r'total |dB/dt| [nT/s]')
        plt.plot(time, np.abs(1e9 * np.sqrt(dB_dt_north_ionosphere**2 + dB_dt_east_ionosphere**2) ), label = r'ionospheric |dB/dt| [nT/s]')
        plt.plot(time, np.abs(1e9 * np.sqrt(dB_dt_north_inner**2 + dB_dt_east_inner**2) ), label = r'inner |dB/dt| [nT/s]')
        plt.plot(time, np.abs(1e9 * np.sqrt(dB_dt_north_outer**2 + dB_dt_east_outer**2) ), label = r'outer |dB/dt| [nT/s]')
        plt.ylim([-8, 8])
        plt.legend()
        filename =  '/wrk-vakka/users/horakons/carrington/plots/{}/GIC/component_dB_dt_timeseries_lat_{}_{}'.format(run,int(lat_deg[i]),run)
        mkdir_path(filename)
        plt.savefig(filename)
        plt.close()
    except:
        print("can't plot components!") 
    '''
    # vv   NEW PLOTS   vv
    # geoelectic field and dB/dt

    fig, (ax1, ax2) = plt.subplots(2, 1)
    plt.rcParams.update({'font.size': 16})

    i0 = 500   # t0 = i0 + 500 seconds for FHA
    #ax1.plot(time[i0:], 1e3 * E_north_arr[ind_min,i0:], label = r'$E_y$ (northward)', color = 'blue')
    #ax1.plot(time[i0:], 1e3 * E_north_arr[ind_min,i0:], label = r'$-E_\theta$', color = 'blue')
    ax1.plot(time[i0:], 1e3 * E_north[i0:], label = r'$E_\lambda$', color = 'blue')
    #ax1.plot(time[i0:], 1e3 * E_north_arr[ind_min,i0:], label = r'$E_\lambda$', color = 'blue')
    #ax1.set_title('Flux Rope {} footpoint, at '.format(int(i+1)) + 
    ax1.set_title('Ionospheric fields, at '.format(int(i+1)) + 
                  r'$\phi$=' + '{:.2f}'.format(phi_deg[i]) + r'$^\circ$, ' + 
                  r'$\lambda$='+ '{:.2f}'.format(lat_deg[i]) + r'$^\circ$' )
    #ax1.set_title('Lat. = {} deg., noon'.format(int(lat_deg[i])))
    #ax1.set_xlabel('time [s]')
    #ax1.set_ylabel(r'$\mathbf{E}$ [V/km]', color = 'blue')
    ax1.set_ylabel(r'$E_\lambda$ [V/km]', color = 'blue')
    ax1.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
    ax1.set_ylim([-0.12, 0.12])
    ax1b = ax1.twinx()
    ax1b.plot(time[i0:], -1e9 * dB_dt_east[i0:], label = r'$-dB_\phi/dt$', color = 'red')
    #ax1b.plot(time[i0:], -1e9 * dB_dt_east[i0:], label = r'$-dB_x/dt$ (westward) [nT/s]', color = 'red')
    ax1b.set_ylabel(r'$-dB_\phi/dt$ [nT/s]', color = 'red')
    ax1b.set_ylim([-1, 1])
    #ax1.legend(loc = 'upper right')
    #ax1b.legend(loc = 'lower right')

    ax2.plot(time[i0:], 1e3 * E_east[i0:], label = r'$E_\phi$', color = 'blue')
    #ax2.plot(time[i0:], 1e3 * E_east_arr[ind_min,i0:], label = r'$E_\phi$', color = 'blue')
    #ax2.plot(time[i0:], 1e3 * E_east_arr[ind_min,i0:], label = r'$E_x$ (eastward)', color = 'blue')
    ax2.set_xlabel('time [s]')
    #ax2.set_ylabel(r'$\mathbf{E}$ [V/km]', color = 'blue')
    ax2.set_ylabel(r'$E_\phi$ [V/km]', color = 'blue')
    ax2.set_yticks([-0.1, -0.05, 0, 0.05, 0.1])
    ax2.set_ylim([-0.12, 0.12])
    ax2b = ax2.twinx()
    ax2b.plot(time[i0:], 1e9 * dB_dt_north[i0:], label = r'$dB_\lambda/dt$', color = 'red')
    #ax2b.plot(time[i0:], 1e9 * dB_dt_north[i0:], label = r'$-dB_\theta/dt$ (northward)', color = 'red')
    #ax2b.plot(time[i0:], 1e9 * dB_dt_north[i0:], label = r'$dB_y/dt$ (northward) [nT/s]', color = 'red')
    #ax2b.set_ylabel(r'$d\mathbf{B}/dt$ [nT/s]', color = 'red')
    ax2b.set_ylabel(r'$dB_\lambda/dt$ [nT/s]', color = 'red')
    ax2b.set_ylim([-1, 1])
    #ax2.legend(loc = 'upper right')
    #ax2b.legend(loc = 'lower right')
    filename =  '/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/FTE_GIC_paper_plots/E_dBdt_lat{}_phi{}_{}'.format(int(lat_deg[i]), int(phi_deg[i]),run)
    mkdir_path(filename)
    plt.savefig(filename)
    plt.close()


