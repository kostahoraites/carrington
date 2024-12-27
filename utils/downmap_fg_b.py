import pytools as pt
import numpy as np
import scipy
from copy import deepcopy

global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations
global R_IONO
R_IONO = R_EARTH + 1e5       # nominal ionosphere altitude 100km (also assumed in Vlasiator)

def theta2lat(theta):
    return -theta + (np.pi / 2)

def lat2theta(lat):
    return (np.pi / 2) - lat

def cartesian_to_spherical(x, y, z):
    '''
    r > 0
    0 < theta < pi
    -pi < phi < pi
    all are assumed to be numpy arrays of equal dimensions

    returns:  r, theta, phi  [tuple]
    '''
    r = (x**2 + y**2 + z**2)**0.5
    theta = np.arccos( z / r )
    phi = np.arctan2(y, x)
    return r, theta, phi

def spherical_to_cartesian(r, theta, phi):
    '''
    r > 0
    0 < theta < pi
    -pi < phi < pi
    all are assumed to be numpy arrays of equal dimensions

    returns:  x, y, z   [tuple]
    '''
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def interpolator_list_3d(f, vec_variable, bounds_error = True, fill_value = None,
                         xmin = None, xmax = None, ymin = None, ymax = None, zmin = None, zmax = None ):
    # vec_variable is a [nx,ny,nz,3] array
    # Create x, y, and z coordinates:
    # Read cellids in order to sort variables
    #cellids = vlsvReader.read_variable("CellID")
    xsize = vec_variable.shape[0]
    ysize = vec_variable.shape[1]
    zsize = vec_variable.shape[2]
    if xmin is None:
        xmin = f.read_parameter('xmin')
    if xmax is None:
        xmax = f.read_parameter('xmax')
    if ymin is None:
        ymin = f.read_parameter('ymin')
    if ymax is None:
        ymax = f.read_parameter('ymax')
    if zmin is None:
        zmin = f.read_parameter('zmin')
    if zmax is None:
        zmax = f.read_parameter('zmax')
    sizes = np.array([xsize, ysize, zsize])
    maxs = np.array([xmax, ymax, zmax])
    mins = np.array([xmin, ymin, zmin])
    #dcell = (maxs - mins)/(sizes.astype('float')).flatten()
    dcell = (maxs - mins)/sizes.astype('float').flatten()
    x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
    y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
    z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
    # Debug:
    if( len(x) != sizes[0] ):
        print("SIZE WRONG: " + str(len(x)) + " " + str(sizes[0]))
    # Create grid interpolation
    interpolator_face_0 = scipy.interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), vec_variable[:,:,:,0], bounds_error = bounds_error, fill_value = fill_value)
    interpolator_face_1 = scipy.interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), vec_variable[:,:,:,1], bounds_error = bounds_error, fill_value = fill_value)
    interpolator_face_2 = scipy.interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), vec_variable[:,:,:,2], bounds_error = bounds_error, fill_value = fill_value)
    return interpolator_face_0, interpolator_face_1, interpolator_face_2

def trace_static_field( f, coord_list, max_iterations, dx, direction='+', fg_b = None, r_trace = 5 * R_EARTH ):
    ''' trace_static_field() integrates along the (static) magnetic field to calculate a final position 
        based on static_field_tracer3d, which in turn is based on Analysator's static_field_tracer

        :param f:         An open vlsv file
        :param coord_list:         List of 1D numpy arrays representing starting positions [x, y, z]
        :param max_iterations:     The maximum amount of iterations before the algorithm stops
        :param dx:                 One iteration step length
        :param direction:          '+' or '-' or '+-' Follow field in the plus direction or minus direction
        :param bvar:               String, variable name to trace [default 'B']
        :returns:                  List of 1D numpy arrays specifying traced positions [x,y,z] 
    '''
    #    if (bvar is not 'B'):
    #        warnings.warn("User defined tracing variable detected. fg, volumetric variable results may not work as intended, use face-values instead.")
    # Read face_B (denoted 'fg_b' in .vlsv files):
    if (fg_b is None):
        fg_b = f.read_variable('fg_b')        # EGL: fg_b.shape = (1024, 736, 736, 3)  
    # Create x, y, and z coordinates:
    # Read cellids in order to sort variables
    interpolators = interpolator_list_3d(f, fg_b, bounds_error = False, fill_value = np.nan)
    if direction == '-':
        multiplier = -1
    else:
        multiplier = 1
    N = coord_list[0].size
    x_out = deepcopy(coord_list[0]); y_out = deepcopy(coord_list[1]); z_out = deepcopy(coord_list[2])
    x_out = x_out.reshape(N)
    y_out = y_out.reshape(N)
    z_out = z_out.reshape(N)          # *_out are 1d arrays
    points = [np.array([x_out, y_out, z_out]).T.reshape([N, 3])]                  # list of 2d arrays each with shape [N, 3]
    point = points[0]
    B_unit = np.zeros([3, N])                                     # B_unit has shape [3,N] for speed considerations
    tf_r_lt = (x_out**2 + y_out**2 + z_out**2)**0.5 > r_trace
    crossed_r_trace = np.zeros(N, dtype = bool)
    for i in range(max_iterations):
        # move along the field line
        B_unit[0, :] = interpolators[0](point)
        B_unit[1, :] = interpolators[1](point)
        B_unit[2, :] = interpolators[2](point)
        B_mag = np.linalg.norm(B_unit, axis=(0))
        B_unit[0, :] = B_unit[0, :] / B_mag
        B_unit[1, :] = B_unit[1, :] / B_mag
        B_unit[2, :] = B_unit[2, :] / B_mag
        new_point = point + multiplier*B_unit.T * dx
        # test if radius r_trace was crossed, update the output
        r_temp = np.linalg.norm(new_point, axis=(1))
        tf_r_lt_new = r_temp > r_trace
        mask = np.logical_and( np.logical_or( np.logical_and( ~tf_r_lt_new, tf_r_lt ), np.logical_and( tf_r_lt_new, ~tf_r_lt ) ),
                              ~crossed_r_trace )      # test for which lines crossed r_trace this iteration, 
                                                      # requiring that r_trace has not been crossed before
        crossed_r_trace[mask] = True
        x_out[mask] = new_point[mask, 0]                # save where the field line crossed radius r_trace
        y_out[mask] = new_point[mask, 1]
        z_out[mask] = new_point[mask, 2]
        tf_r_lt = tf_r_lt_new
        point = new_point
    x_out[~crossed_r_trace] = point[~crossed_r_trace, 0]     # for remaining points that never cross the threshold, save their final location
    y_out[~crossed_r_trace] = point[~crossed_r_trace, 1]
    z_out[~crossed_r_trace] = point[~crossed_r_trace, 2]
    x_out = x_out.reshape(coord_list[0].shape)
    y_out = y_out.reshape(coord_list[0].shape)
    z_out = z_out.reshape(coord_list[0].shape)
    crossed_r_trace = crossed_r_trace.reshape(coord_list[0].shape)
    #    return points 
    return x_out, y_out, z_out, crossed_r_trace

def trace_coord(f, coord_list, fg_b = None, trace_method='dipole', direction = '+', coord_in='spherical', coord_out='cartesian', r_trace = 5 * R_EARTH, max_iterations = 1000, dx = R_EARTH / 50):
    '''
    trace_coord() takes input coordinates given by field_list and maps them to a location in the magnetosphere by following trajectories along field lines.
    
    f: a VlsvReader object OR the .vlsv file name.
    coord_list: a 3-element list of numpy arrays of equal length (or scalars?),
                specifying the [SI] position to start tracing from
                ex. [ [1.2e7 2.3e7 3.1e7] , [0 0 0] , [1e7 2e7 3e7] ]
    fg_b: if set to None, then fg_b is read from f in the program
          otherwise it should be set to the output of f.read_fsgrid_variable('fg_b') when the function is called
    trace_method: the method used to follow the particles along the field lines ---  'dipole', 'integrateB', 'particlePusher'
    coord_in: 'spherical' or 'cartesian' the coordinate system (r, theta, phi) or (x, y, z) [GSE]
    coord_out: 'spherical' or 'cartesian', the coordinate system of the output (locations where the field lin)
    r_trace: The radius to trace the field lines to.
             In practice, the minimum radius of the vlasov grid [meters] (vg_rmin) default value is found by hand in EGL run.
    Note: spherical means r>0, 0<theta<pi, (phi doesn't matter but assume -pi < phi < pi)
 
    '''
    if type(f) == str: #convert filename to .vlsv object
        f = pt.vlsvfile.VlsvReader(f)
    if (fg_b is None) and (trace_method != 'dipole'):
        fg_b = f.read_fsgrid_variable('fg_b')        # EGL: fg_b.shape = (1024, 736, 736, 3)  
    if coord_in == 'spherical':
        r, theta, phi = coord_list[0], coord_list[1], coord_list[2]
        x, y, z = spherical_to_cartesian(r, theta, phi)
    elif coord_in == 'cartesian':
        x, y, z =  coord_list[0], coord_list[1], coord_list[2]
        r, theta, phi = cartesian_to_spherical(x, y, z)
    if trace_method == 'dipole':
        # r = L cos^2 (lat), where L is the "L-shell" (magnetic dipole)
        # keeping the L-shell constant (moving along a field line): r1 / r2 = cos^2(lat1) / cos^2(lat2)
        # find the latitude where field lines cross r_trace
        lat1 = theta2lat(theta)                      # convert theta to latitude (radians)
        lat2 = ( abs(lat1) / lat1) * np.arccos( np.cos(lat1) * (r_trace / r)**0.5 )            # first term on RHS accounts for the sign
        theta_out = lat2theta(lat2)
        x_out, y_out, z_out = spherical_to_cartesian( x*0 + r_trace, theta_out, phi)           # phi is assumed to stay the same    
    elif trace_method == 'integrateB':
        x_out, y_out, z_out, ind = trace_static_field( f, [x,y,z], max_iterations, dx, direction=direction, fg_b = fg_b, r_trace = r_trace )
        x_out[~ind] = np.nan
        y_out[~ind] = np.nan
        z_out[~ind] = np.nan
    if coord_out == 'spherical':
        r_out, theta_out, phi_out = cartesian_to_spherical(x_out, y_out, z_out)
        coord_list_out = [r_out, theta_out, phi_out]
    elif coord_out == 'cartesian':
        coord_list_out = [x_out, y_out, z_out]
    return coord_list_out



# Downmap to the ionosphere
filename = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk_with_fg_10/bulk_with_fg_10.0000101.vlsv'

R_INNER = 4.7 * R_EARTH  # inner boundary
dx = 1e4  # meters 
max_iterations = int(30 * R_EARTH / dx)

f = pt.vlsvfile.VlsvReader(filename)
fg_b = f.read_variable('fg_b')

#Input points (magnetosphere):
coord_list_msphere = np.array([[-10*R_EARTH, -11*R_EARTH], [0,0], [1*R_EARTH,1*R_EARTH]])  # [x_list, y_list, z_list]

#trace magnetosphere --> R_INNER (tracing fg_b fields):
coord_list_traced_tmp = trace_coord(filename, coord_list_msphere, fg_b = fg_b, trace_method = 'integrateB', direction = '+',
                                    coord_in = 'cartesian', coord_out = 'cartesian', r_trace = R_INNER, max_iterations = max_iterations, dx = dx)
#trace R_INNER --> R_IONO (dipole formula):
coord_list_traced_final = trace_coord(filename, coord_list_traced_tmp, fg_b = fg_b, trace_method = 'dipole', direction = '+',
                                    coord_in = 'cartesian', coord_out = 'cartesian', r_trace = R_IONO, max_iterations = max_iterations, dx = dx)

print('starting coordinates:\n', np.array(coord_list_msphere))
print('downmapped coordinates:\n', np.array(coord_list_traced_final))
