# trace the magnetic field that connects to the 

import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import matplotlib.colors as colors
import numpy as np
import pytools as pt
from myutils import cartesian_to_spherical, spherical_to_cartesian, mkdir_path, timer, timeseries, get_vlsvfile_fullpath
import scipy
import statsmodels.api as sm
from fieldtracer import static_field_tracer_3d
import argparse
from matplotlib import ticker

parser = argparse.ArgumentParser()
parser.add_argument('-nproc', default=1, help="number of processors to use " )
global ARGS
ARGS = parser.parse_args()

global f0
global R_E
R_E = 6371000.
#run = "FHAFGB"   # FHA files that contain full resolution 'fg_b' magnetic field variable
run = "FHA"   # FHA files that contain full resolution 'fg_b' magnetic field variable
fileIndex = 1100 # 900, 1000, 1100, 1165  # 110  # time = fileIndex*10 for FHAFGB run (larger files saved more sparsely).
f0 = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(run, fileIndex))   # time in Fig. 1: t = 1165 s
#pos = f0.get_ionosphere_node_coords()  # shape (21568, 3)

#f_iono = pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/ig_B/ionosphere_B_sidecar_FHA.0000784.vlsv")
#pos = f_iono.read_variable('ig_r')   # ionospheric grid.  array dimensions (43132, 3)

# Load grid coordinates

#vg_coordinates = f0.read_variable("vg_coordinates")
#vg_r, vg_theta, vg_phi = cartesian_to_spherical(vg_coordinates[:,0], vg_coordinates[:,1], vg_coordinates[:,2])
ig_coordinates = f0.get_ionosphere_node_coords()  # node locations (element corners): shape (21568, 3)
ig_r, ig_theta, ig_phi = cartesian_to_spherical(ig_coordinates[:,0], ig_coordinates[:,1], ig_coordinates[:,2])
ig_theta_deg = ig_theta * 180. / np.pi
ig_phi_deg = ig_phi * 180. / np.pi

# Coordinate of interest (ig_)
'''
ind_ig_plot = 8       # [0.12, -0.05, 1.01] R_E. This node shows quasiperiodic ig_fac, with 3 clear spikes between 1100-1598s (FHA run)
ig_testcoord = ig_coordinates[ind_ig_plot, :]
'''

# trace field lines starting from ionospheric coordinates theta, phi
#dayside aurora

ind = np.where( (np.abs(ig_phi_deg) < 30) & 
                (np.abs(90 - ig_theta_deg) > 77) &
                (np.abs(90 - ig_theta_deg) < 83) )[0]

#for i in ind:
def trace_it(i):
    # Find nearest neighbor of the ionosphere grid, index by 'ind', to the specified lat and phi
    save = False

    print('index ', i, ', ', ind.size, ' total')

    # Initial coordinate in the vg grid (upmapped from ig_testcoord)
    ig_upmappednodecoords = f0.read_variable('ig_upmappednodecoords')
    vg_coordinates_0 = ig_upmappednodecoords[i, :]

    filename = '/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/txt_files/btraces/fg_innerboundary_btrace_C2_{}_t{}_phi{:.2f}_theta{:.2f}.txt'.format(run, int(fileIndex), ig_phi_deg[i], ig_theta_deg[i])
    if save:
        dx = 1e4
        # Trace the B-field line
        if ig_theta_deg[i] > 90:
            direction = '+'
        else:
            direction = '-'
        # Set parameters of the tracing and output file
        max_iterations = 6000
        #ncoordsave = 150    # total number of coordinates in the final file
        #dstep_write = int(max_iterations / ncoordsave)   # number of iterations between writes  (*2 if keyword direction = '+-')
        x_b = static_field_tracer_3d( f0, np.array([vg_coordinates_0]), max_iterations, dx, direction=direction, grid_var='vg_b_vol' )   # numpy array
        # Save to a text file
        x = x_b[0, :,:]
        np.savetxt(filename, x)
    else:
        # now load the data and save a variable for plotting:
        x=np.loadtxt(filename)

    dx = np.linalg.norm(x[1,:] - x[0,:])

    npts = x.shape[0]  # number of traced points

    tmin = 501
    tmax = 1612  # 510
    nt = tmax - tmin + 1

    x = x[:, :]

    xvar = np.array(nt * [np.arange(npts) * dx] ).transpose() / R_E
    tvar = np.arange(tmin, tmax+1)[None, :] * (np.zeros([npts, nt])+1)

    #LOAD DATA

    x_0_str = 't_{}_ig_{}_x{:.2f}_y{:.2f}_z{:.2f}_RE'.format(int(fileIndex), i, x[0,0]/R_E, x[0, 1]/R_E, x[0, 2]/R_E)
    #x_0_str = 'ig_{}_t{}_x{:.2f}_y{:.2f}_z{:.2f}_RE'.format(i, int(fileIndex), x[0,0]/R_E, x[0, 1]/R_E, x[0, 2]/R_E)
    datafile = 'txtfiles/J_par_B'+ x_0_str +'.txt'
    datafile_xo = 'txtfiles/XO'+ x_0_str +'.txt'
    if save:
        zvar = np.zeros([npts, nt])
        xo_var = np.zeros([npts, nt])
        for t_i in range(tmin, tmax+1):
            f = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(run, t_i))
            f_xo = pt.vlsvfile.VlsvReader('/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/pysidecar_sdf_bulk1.{}.vlsv'.format(str(t_i).zfill(7)))
            B = f.read_interpolated_variable('vg_b_vol', x)
            J = f.read_interpolated_variable('vg_j', x)
            B_mag = np.linalg.norm(B, axis = 1)
            J_par = np.sum(J * B, axis = 1) / B_mag
            zvar[:, t_i - tmin] = J_par / B_mag   # units: H^-1
            try:
                vg_lmn_nld = f_xo.read_interpolated_variable('vg_lmn_neutral_line_distance', x, method = 'nearest')   # <0.2 to find XO lines (what are the units?)
                vg_dBNdL = f_xo.read_interpolated_variable('vg_dBNdL', x, method = 'nearest')                         # -3e-16 to 3e-16 
                xo_var[:, t_i - tmin] = (vg_lmn_nld < 0.2) * np.sign(vg_dBNdL)              # -1: 0-line, 1: X-line, 0: neither
            except:
                print('Data missing! time step= {}'.format(t_i))
        np.savetxt(datafile, zvar)
        np.savetxt(datafile_xo, xo_var, fmt='%d')
    else:
        zvar = np.loadtxt(datafile)
        xo_var = np.loadtxt(datafile_xo)

    # prune the data
    t0 = 900 # 1000 # starting plot time
    tf = 1612
    i_t0 = t0 - tmin
    i_tf = tf - tmin
    i_x0 = 0 # starting plot position index
    i_xf = 5000 + 1
    npts_raw_data = i_xf - i_x0
    # Alfven wave x(t) propagation
    va_1d = f0.read_interpolated_variable('vg_va', x)
    x_1d = np.arange(npts_raw_data) * dx
    t_1d = np.cumsum(dx / va_1d)
    t_1d_reverse = np.cumsum(dx / np.flip(va_1d))

    dx_cell = 1e6
    di_cell = int(dx_cell/dx)
    print('di_cell', di_cell)
    x = x[i_x0:i_xf:di_cell, :]
    xvar = xvar[i_x0:i_xf:di_cell, i_t0:i_tf]
    tvar = tvar[i_x0:i_xf:di_cell, i_t0:i_tf]
    zvar = zvar[i_x0:i_xf:di_cell, i_t0:i_tf]
    xo_var = xo_var[i_x0:i_xf:di_cell, i_t0:i_tf]
    x_1d = x_1d[i_x0:i_xf:di_cell]
    t_1d = t_1d[i_x0:i_xf:di_cell]
    t_1d_reverse = t_1d_reverse[i_x0:i_xf:di_cell]
    npts_plot = x.shape[0]


    # O-line occurrences
    i_Os, j_Os = np.where(xo_var == -1)

    #PLOT
    plt.rcParams.update({'font.size': 16})
    cmap = 'RdBu'  # 'bwr', 'plasma'

    low, high = (np.nanmin(zvar), np.nanmax(zvar))
    # find the new limits
    bound = max(abs(low), abs(high))
    #if bound > 1:
    #    bound = 0.01  # default

    fig, (ax, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(12, 6.5)
    #linthresh = bound /1000.

    # First plot: FAC keogram
    linthresh = np.max(np.abs(zvar[0,:]))*2
    im = ax.pcolormesh(xvar, tvar, zvar, cmap = cmap, shading = 'nearest', 
                       norm=colors.SymLogNorm(linthresh=linthresh, linscale=1, vmin=-bound, vmax=bound, base=10) ) # norm = norm
    ax.plot(x_1d/R_E, fileIndex + t_1d, color = 'yellow', linewidth = 3., linestyle = '--')
    ax.plot(np.flip(x_1d)/R_E, fileIndex + t_1d_reverse, color = 'lime', linewidth = 3., linestyle = '--')
    ax.set_xlabel(r'dist. along curve [$R_E$]')
    ax.set_ylabel('time [s]')
    ax.set_title('FACs on cusp field lines')
    #ax.set_xlim([0, npts_raw_data*dx])

    # set new axes and colorbar limits
    ax.set_ylim([t0, tf])
    im.set_clim([-bound, bound])
    #annotate endpoints
    ax.annotate(r'$x_3$',
                xy=(0,t0-(tf-t0)/12.), xytext=(0,t0-(tf-t0)/12.),
                annotation_clip=False, rotation = 0, color = 'orange', fontsize = 26)
    #ax.annotate('[{:.2f},{:.2f},{:.2f}]'.format(x[0, 0]/R_E, x[0, 1]/R_E, x[0, 2]/R_E) + r' $R_E$',
    #            xy=(0,t0-(tf-t0)/2.), xytext=(0,t0-(tf-t0)/2.),
    #            annotation_clip=False, rotation = -30., color = 'orange')
    ax.annotate(r'$x_4$',
                xy=(npts_raw_data*dx/R_E,t0-(tf-t0)/12.), xytext=(npts_raw_data*dx/R_E,t0-(tf-t0)/12.),
                annotation_clip=False, rotation = 0, color = 'orange', fontsize = 26)
    #ax.annotate('[{:.2f},{:.2f},{:.2f}]'.format(x[-1, 0]/R_E, x[-1, 1]/R_E, x[-1, 2]/R_E) + r' $R_E$', 
    #            xy=(npts_raw_data*dx/R_E,t0-(tf-t0)/2.), xytext=(npts_raw_data*dx/R_E,t0-(tf-t0)/2.),
    #            annotation_clip=False, rotation = -30., color = 'orange')
    #annotate O-lines
    for i_O, j_O in zip(i_Os, j_Os):
        ax.annotate('O', xy=(xvar[i_O, j_O], tvar[i_O, j_O]), xytext=(xvar[i_O, j_O], tvar[i_O, j_O]), annotation_clip=False)

    #ax.set_position([1,1,4,5], which = 'both')

    #color bar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)

    # tick marks

    maxlog=int(np.floor(np.log10(bound)))
    minlog = maxlog
    logthresh=int(np.floor(np.log10(linthresh)))
    logstep=1
    def signstr(number):
        if number > 0:
            return '^+'
        else:
            return '^-'
    ticks=([-(10**xx) for xx in range(logthresh, minlog+1, logstep)][::-1]
            +[0.0]
            +[(10**xx) for xx in range(logthresh, maxlog+1, logstep)] )
    ticklabels=([r'$-10^{}$'.format('{' + str(xx) + '}') for xx in range(logthresh, minlog+1, logstep)][::-1]
            +['0']
            +[r'$10^{}$'.format('{' + str(xx) + '}') for xx in range(logthresh, maxlog+1, logstep)])
    allticks=([-(10**(xx/2.)) for xx in range(logthresh-1, 2*minlog+1, logstep)][::-1]
            +[0.0]
            +[(10**(xx/2.)) for xx in range(logthresh-1, 2*maxlog+1, logstep)] )
    allticks_ind = np.where(np.abs(np.array(allticks)) < bound)[0]
    allticks = np.array(allticks)[allticks_ind]
    print('ticks: ', ticks)
    #cbformat = ticker.ScalarFormatter()
    #cbformat.set_scientific('%.2e')
    #cbformat.set_powerlimits((-10,10))
    #cbformat.set_useMathText(True)
    mycbar = fig.colorbar(im, cax=cax, orientation='vertical') #, ticks = ticks) #, format = cbformat)
    mycbar.set_label(r'$J_\parallel/B$ $[H^{-1}]$')
    #mycbar.set_ticklabels(ticklabels)
    mycbar.ax.minorticks_on()
    mycbar.ax.yaxis.set_ticks(allticks, minor=True)
    #plt.tight_layout()

    # Second plot: magnetosphere cut

    bulkname = get_vlsvfile_fullpath(run, fileIndex)
    r_C = 5.6 # coupling radius [R_E]. calculated from f.get_config()['ionosphere']['downmapRadius'][0]

    cmap = matplotlib.cm.plasma
    cmap.set_under('grey')

    outputdir = '/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/FTE_GIC_paper_plots/'
    pt.plot.plot_colormap3dslice(filename=bulkname,var='proton/vg_pressure', boxre=[0, 14, -14, 14], normal = 'y', run=run,
                            colormap='plasma',step=fileIndex,outputdir=outputdir, vmin=0,vmax=1e-9,lin = True,
                            Earth=1, streamlines='vg_b_vol', streamlinedensity=2.5, streamlinethick = 2, streamlinecolor = 'black',
                            cutpointre=x[-1,1,]/R_E, axes = ax2, scale=2.2, useimshow=True, draw = 1)
    ax2.plot(x[:,0]/R_E, x[:,2]/R_E, color = 'lime', linewidth=5, alpha = 0.7)
    ax2.scatter(x[0,0]/R_E, x[0,2]/R_E, color = 'lime', s= 200, alpha = 0.7)
    ax2.scatter(x[-1,0]/R_E, x[-1,2]/R_E, color = 'lime', s = 200, alpha = 0.7)
    ax2.annotate(r'$x_3$', xy=(x[0,0]/R_E, x[0,2]/R_E), xytext=(x[0,0]/R_E, x[0,2]/R_E), color = 'white', annotation_clip=False, fontsize = 20)
    ax2.annotate(r'$x_4$', xy=(x[-1,0]/R_E, x[-1,2]/R_E), xytext=(x[-1,0]/R_E, x[-1,2]/R_E), color = 'white', annotation_clip=False, fontsize = 20)
    #ax2.set_position([5.5,1,4,5], which = 'both')

    theta = np.arange(0, 2*np.pi, 0.01)
    ax2.plot(r_C * np.cos(theta), r_C * np.sin(theta), linestyle = '--', color = 'grey', zorder = 2, linewidth=3)

    # save figure
    filename = '/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/FTE_GIC_paper_plots/keograms/{}/FACs_Jpar_keogram'.format(int(fileIndex))+'_ig_{}_'.format(i)+x_0_str+'.png'
    mkdir_path(filename)
    plt.tight_layout()
    plt.savefig(filename, bbox_inches="tight")
    plt.close()


## v   Parallel processing     v

#UNCOMMENT TO RECALCULATE TSYGANENKO OCB

from multiprocessing import Pool
pool = Pool(int(ARGS.nproc))
ocb_control = pool.map(trace_it, ind)
pool.close()
pool.join()



