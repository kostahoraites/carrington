# trace the magnetic field that connects to the 

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pytools as pt
#from numba import jit
from myutils import cartesian_to_spherical, spherical_to_cartesian, mkdir_path, timer, timeseries, get_vlsvfile_fullpath
import scipy
import statsmodels.api as sm
from fieldtracer import static_field_tracer_3d

global R_EARTH
R_EARTH = 6371000.
run = "FHAFGB"   # FHA files that contain full resolution 'fg_b' magnetic field variable
fileIndex = 110  # time = fileIndex*10 for FHAFGB run (larger files saved more sparsely).

f = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(run, fileIndex))

# Load grid coordinates
vg_coordinates = f.read_variable("vg_coordinates")
vg_r, vg_theta, vg_phi = cartesian_to_spherical(vg_coordinates[:,0], vg_coordinates[:,1], vg_coordinates[:,2])
#ig_coordinates = dct_E_sidecar['ig_r']   # this is locations of triangular elements (barycenters), NOT the same as node locations
ig_coordinates = f.get_ionosphere_node_coords()  # node locations (element corners): shape (21568, 3)

# Coordinate of interest (ig_)
ind_ig_plot = 8       # [0.12, -0.05, 1.01] R_E. This node shows quasiperiodic ig_fac, with 3 clear spikes between 1100-1598s (FHA run)
ig_testcoord = ig_coordinates[ind_ig_plot, :]

# Initial coordinate in the vg grid (upmapped from ig_testcoord)
ig_upmappednodecoords = f.read_variable('ig_upmappednodecoords')
vg_coordinates_0 = ig_upmappednodecoords[ind_ig_plot, :]

# Set parameters of the tracing and output file
max_iterations = 6000
dx = 1e4
ncoordsave = 150    # total number of coordinates in the final file
dstep_write = int(max_iterations / ncoordsave)   # number of iterations between writes  (*2 if keyword direction = '+-')

# Trace the B-field line
x_b = static_field_tracer_3d( f, [vg_coordinates_0], max_iterations, dx, direction='-', fg = 'fg_b' )   # list of 3-element lists

# Save to a text file
output_file = open('txt_files/fg_innerboundary_btrace_{}.txt'.format(run), 'w')
i = 0
for x in x_b:
    if i == dstep_write:
        output_file.write(str(x[0]/R_EARTH).replace('[','').replace(']','') + '\n')
        i=0
    else:
        i+=1

output_file.close()


