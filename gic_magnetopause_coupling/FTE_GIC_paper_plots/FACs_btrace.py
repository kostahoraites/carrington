# trace the magnetic field that connects to the 

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytools as pt
from myutils import cartesian_to_spherical, spherical_to_cartesian, mkdir_path, timer, timeseries, get_vlsvfile_fullpath
import scipy
import statsmodels.api as sm
#from static_field_tracer_3d_alt import static_field_tracer_3d_alt
from fieldtracer import static_field_tracer_3d

global R_EARTH
R_EARTH = 6371000.
#run = "FHAFGB"   # FHA files that contain full resolution 'fg_b' magnetic field variable
run = "FHA"   # FHA files that contain full resolution 'fg_b' magnetic field variable
fileIndex = 1165  # 110  # time = fileIndex*10 for FHAFGB run (larger files saved more sparsely).
f = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(run, fileIndex))   # time in Fig. 1: t = 1165 s
pos = f.get_ionosphere_node_coords()  # shape (21568, 3)

#f_iono = ft.f("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/ig_B/ionosphere_B_sidecar_FHA.0000784.vlsv")
#pos = f_iono.read_variable('ig_r')   # ionospheric grid.  array dimensions (43132, 3)



# Load grid coordinates

vg_coordinates = f.read_variable("vg_coordinates")
vg_r, vg_theta, vg_phi = cartesian_to_spherical(vg_coordinates[:,0], vg_coordinates[:,1], vg_coordinates[:,2])
ig_coordinates = f.get_ionosphere_node_coords()  # node locations (element corners): shape (21568, 3)

# Coordinate of interest (ig_)
'''
ind_ig_plot = 8       # [0.12, -0.05, 1.01] R_E. This node shows quasiperiodic ig_fac, with 3 clear spikes between 1100-1598s (FHA run)
ig_testcoord = ig_coordinates[ind_ig_plot, :]
'''

#f = ft.f(get_vlsvfile_fullpath('FHA', 1165))

# calculate geoelectric field at footpoints of Figure 1 flux ropes

#magnetospheric coord1, coord2

coord1 = np.array([10.450000, -1.430000, -2.600000]) * R_EARTH    # 'FLUX ROPE 1'
coord2 = np.array([10.500000, 1.470000, -2.560000]) * R_EARTH     # 'FLUX ROPE 2'

#calculate ionospheric downmapped coordinates

xyz1 = f.read_interpolated_variable('vg_connection_coordinates_bw', coord1)/R_EARTH
xyz2 = f.read_interpolated_variable('vg_connection_coordinates_bw', coord2)/R_EARTH
r1, theta1, phi1 = cartesian_to_spherical( *xyz1)
r2, theta2, phi2 = cartesian_to_spherical( *xyz2)

# select just the field line rooted at coord2 ('Flux rope 2' in the paper)

theta = theta2  # np.array([theta1, theta2])
phi = phi2     # np.array([phi1, phi2])
'''
lat_deg = np.array([-79., -80., -81., 79., 80., 81.])
lat = lat_deg * np.pi / 180.   # (auroral) latitude
theta = (np.pi / 2) - lat  # co-latitude
phi = lat * 0 + 0.00000000000001              # noon, klug to make phi positive
'''

#theta_deg = theta * 180. / np.pi
#lat_deg = 90. - theta_deg
#phi_deg = phi * 180. / np.pi

x0, y0, z0 = spherical_to_cartesian(R_EARTH, theta, phi)



# Find nearest neighbor of the ionosphere grid, index by 'ind_ig_plot', to the specified lat and phi
dist = np.sqrt((x0 - pos[:,0])**2 + (y0 - pos[:,1])**2 + (z0 - pos[:,2])**2)
ind_ig_plot = np.argmin(dist)



# Initial coordinate in the vg grid (upmapped from ig_testcoord)
ig_upmappednodecoords = f.read_variable('ig_upmappednodecoords')
vg_coordinates_0 = ig_upmappednodecoords[ind_ig_plot, :]

# Set parameters of the tracing and output file
max_iterations = 6000
dx = 1e4
ncoordsave = 150    # total number of coordinates in the final file
dstep_write = int(max_iterations / ncoordsave)   # number of iterations between writes  (*2 if keyword direction = '+-')

# Trace the B-field line

x_b = static_field_tracer_3d( f, np.array([vg_coordinates_0]), max_iterations, dx, direction='+', grid_var='vg_b_vol' )   # list of 3-element lists
#x_b = static_field_tracer_3d_alt( f, [vg_coordinates_0], max_iterations, dx, direction='+', fg_b = None )   # list of 3-element lists
print(x_b)
print(type(x_b))

# Save to a text file
filename = '/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/txt_files/fg_innerboundary_btrace_C2_{}.txt'.format(run)
np.savetxt(filename, x_b[0, :,:])

'''
output_file = open(filename, 'w')

i = 0
for x in x_b:
    if i == dstep_write:
        output_file.write(str(x[0]/R_EARTH).replace('[','').replace(']','') + '\n')
        i=0
    else:
        i+=1
'''
output_file.close()


