import ftest as ft
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from myutils import get_vlsvfile_fullpath, mkdir_path

R_E = 6.371e6
#x=np.loadtxt('/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/txt_files/fg_innerboundary_btrace_C2_FHA.txt')
x = R_E * np.loadtxt('/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/virtual_sc_jonas_keogram/fg_innerboundary_btrace_FHAFGB_new.txt')


dx = np.linalg.norm(x[1,:] - x[0,:])

npts = x.shape[0]  # number of traced points


tmin = 501
tmax = 1612  # 510
nt = tmax - tmin + 1

xvar = np.array(nt * [np.arange(npts) * dx] ).transpose() / R_E
tvar = np.arange(tmin, tmax+1)[None, :] * (np.zeros([npts, nt])+1)

run = 'FHA'

#LOAD DATA

f = ft.f(get_vlsvfile_fullpath(run, 1165))
va_1d = f.read_interpolated_variable('vg_va', x)
x_1d = np.arange(npts) * dx
t_1d = np.cumsum(dx / va_1d)

x_0_str = '_x{:.2f}_y{:.2f}_z{:.2f}_RE'.format(x[0,0]/R_E, x[0, 1]/R_E, x[0, 2]/R_E)
save = False
datafile = 'J_par_B'+ x_0_str +'.txt'
if save:
    zvar = np.zeros([npts, nt])
    for t_i in range(tmin, tmax+1):
        f = ft.f(get_vlsvfile_fullpath(run, t_i))
        B = f.read_interpolated_variable('vg_b_vol', x)
        J = f.read_interpolated_variable('vg_j', x)
        B_mag = np.linalg.norm(B, axis = 1)
        J_par = np.sum(J * B, axis = 1) / B_mag
        zvar[:, t_i - tmin] = J_par / B_mag   # units: H^-1
        np.savetxt(datafile, zvar)
else:
    zvar = np.loadtxt(datafile)



#PLOT
plt.rcParams.update({'font.size': 14})
cmap = 'bwr'  # 'bwr', 'plasma'

fig, ax = plt.subplots()
im = ax.pcolormesh(xvar, tvar, zvar, cmap = cmap) # norm = norm
ax.plot(x_1d/R_E, tmin + t_1d, color = 'yellow', linewidth = 3., linestyle = '--')
ax.set_xlabel(r'dist. along curve [$R_E$]')
ax.set_ylabel('time [s]')
ax.set_title('FACs on cusp field lines')
#ax.set_xlim([0, npts*dx])
ax.set_ylim([tmin, tmax])
ax.annotate('[{:.2f},{:.2f},{:.2f}]'.format(x[0, 0]/R_E, x[0, 1]/R_E, x[0, 2]/R_E) + r' $R_E$', xy=(0,tmin-(tmax-tmin)/2.), xytext=(0,tmin-(tmax-tmin)/2.),
            annotation_clip=False, rotation = -30., color = 'orange')
ax.annotate('[{:.2f},{:.2f},{:.2f}]'.format(x[-1, 0]/R_E, x[-1, 1]/R_E, x[-1, 2]/R_E) + r' $R_E$', xy=(npts*dx/R_E,tmin-(tmax-tmin)/2.), xytext=(npts*dx/R_E,tmin-(tmax-tmin)/2.),
            annotation_clip=False, rotation = -30., color = 'orange')
   #color bar
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
mycbar = fig.colorbar(im, cax=cax, orientation='vertical')
mycbar.set_label(r'$J_\parallel/B$ $[H^{-1}]$')
plt.tight_layout()

filename = 'FACs_Jpar_keogram'+x_0_str+'.png'
mkdir_path(filename)
plt.savefig(filename)
plt.close()


