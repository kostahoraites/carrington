#from pyCalculations import fieldtracer

import fieldtracertest
from static_field_tracer_3d_alt import static_field_tracer_3d_alt
from static_field_tracer3d import static_field_tracer3d
from static_field_tracer_3d import static_field_tracer_3d   # this version will be uploaded to analysator
import pytools as pt
import matplotlib.pyplot as plt
import numpy as np

# main program

R_EARTH = 6.371e6            #check what value is used in simulations

f0=pt.vlsvfile.VlsvReader("/wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.0001750.vlsv")
f=pt.vlsvfile.VlsvReader("/wrk/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.0001760.vlsv")
fg_b = f.read_fsgrid_variable('fg_b')        # fg_b.shape = (1024, 736, 736, 3)
dx = R_EARTH / 100

x0 = np.array([R_EARTH / 2**0.5, 0, R_EARTH / 2**0.5])
x_arr = np.array([x0[0], x0[0]])
y_arr = np.array([x0[1], x0[1]])
z_arr = np.array([x0[2], x0[2]])
coord_list = [ x_arr, y_arr, z_arr ]
max_iterations = 4
direction = '-'

points_3d_alt = static_field_tracer_3d_alt( f, coord_list, max_iterations, dx, direction=direction, fg_b = fg_b )
print('points_3d_alt', points_3d_alt)

# return points    # a list of numpy arrays

points_3d = static_field_tracer3d( f, x0, max_iterations, dx, direction=direction, bvar=fg_b )
print('points_3d', points_3d)

points_3d_pm = static_field_tracer3d( f, x0, max_iterations, dx, direction='+-', bvar=fg_b )
print('points_3d_pm', points_3d_pm)

points_test = fieldtracertest.static_field_tracer( f, x0, max_iterations, dx, direction='+-', bvar=fg_b )
print('points_test', points_test)

# return points

vlsvReader_list = [f0, f]
points_dynamic = fieldtracertest.dynamic_field_tracer( vlsvReader_list, x0, max_iterations, dx)

#fieldtracer.static_field_tracer
# static_field_tracer( vlsvReader, x0, max_iterations, dx, direction='+', bvar='B' )
# return points






