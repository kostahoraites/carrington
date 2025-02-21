# CHECK Fukushima's theorem (the geoelectric field should be approximately divergence-free) 

import numpy as np
import scipy
import ftest as ft
from myutils import cartesian_to_spherical
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable



def symlims(data):
    amplitude = np.nanmax((np.abs(data)))
    return (-amplitude, amplitude)

# Load ionosphere E-field data:

fileIndex = 1300
filename = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/geoelectric_field/ionosphere_gic_FHA_{}.vlsv'.format(str(fileIndex).zfill(7))

f = ft.f(filename)

ig_coordinates = f.read_variable('ig_r')
ig_E_theta = -f.read_variable('ig_E_north')
ig_E_phi = f.read_variable('ig_E_east')

ig_r, ig_theta, ig_phi = cartesian_to_spherical(ig_coordinates[:,0], ig_coordinates[:,1], ig_coordinates[:,2])
ig_theta_deg = ig_theta * 180. / np.pi
ig_phi_deg = ig_phi * 180. / np.pi

# Construct interpolating grid:

nphi = 360j
ntheta = 180j

grid_phi, grid_theta = np.mgrid[-np.pi:np.pi:nphi, 0:np.pi:ntheta]
grid_phi_deg = grid_phi * 180/np.pi; grid_theta_deg = grid_theta * 180/np.pi
grid_lat_deg = 90. - grid_theta_deg

# Interpolate:

points = (ig_phi, ig_theta) 
xi = (grid_phi, grid_theta)

E_phi_interp = scipy.interpolate.griddata(points, ig_E_phi, xi, method='cubic')
E_theta_interp = scipy.interpolate.griddata(points, ig_E_theta, xi, method='cubic')

# Plot E-field:

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
lims = symlims(E_theta_interp)
im = ax.pcolormesh(grid_phi_deg, grid_theta_deg, E_theta_interp,
                   cmap='bwr', vmin = lims[0], vmax = lims[1], shading='auto')
ax.set_xlabel(r'$\phi$ [deg.]')
ax.set_ylabel(r'$\theta$ [deg.]')
ax.set_title(r'$E_\theta$')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
mycbar = fig.colorbar(im, cax=cax, orientation='vertical')
mycbar.set_label(r'$V/m$')
ax.streamplot(grid_phi_deg[:,0], grid_theta_deg[0,:], E_phi_interp.transpose(), E_theta_interp.transpose(), color='black', linewidth=1, cmap=plt.cm.inferno, density=4, arrowstyle='->', arrowsize=1.5)
plt.savefig('E_theta_interp.png')
plt.close()


fig, ax = plt.subplots()
lims = symlims(E_phi_interp)
im = ax.pcolormesh(grid_phi_deg, grid_theta_deg, E_phi_interp,
                   cmap='bwr', vmin = lims[0], vmax = lims[1], shading='auto') 
ax.set_xlabel(r'$\phi$ [deg.]')
ax.set_ylabel(r'$\theta$ [deg.]')
ax.set_title(r'$E_\phi$')
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
mycbar = fig.colorbar(im, cax=cax, orientation='vertical')
mycbar.set_label(r'$V/m$')
ax.streamplot(grid_phi_deg[:,0], grid_theta_deg[0,:], E_phi_interp.transpose(), E_theta_interp.transpose(), color='black', linewidth=1, cmap=plt.cm.inferno, density=4, arrowstyle='->', arrowsize=1.5)
plt.savefig('E_phi_interp.png')
plt.close()

# Calculate Curl(E), Div(E):



# Plot Curl(E), Div(E)




