import pytools as pt
from pyCalculations.intpol_points import vlsv_intpol_points
import pyCalculations.fieldtracer as fieldtracer
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from static_field_tracer_3d_alt import static_field_tracer_3d_alt
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path, cartesian_to_spherical, spherical_to_cartesian, numcurl3d, numjacobian3d
import os, sys
import argparse
import warnings
from copy import deepcopy
import scipy
from scipy.optimize import curve_fit
import pandas
from time import time
#import magnetopause3dk          # don't want to do this at the top of the file because it kind of breaks vorna
#from magnetopause3dk import make_streamlines, get_magnetopause


global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations
#global CELLSIZE
global ROOT_DIR
ROOT_DIR = '/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/new_plots/'
global mu_0
mu_0 = 4e-7 * np.pi


def theta2lat(theta):
    return -theta + (np.pi / 2)


def lat2theta(lat):
    return (np.pi / 2) - lat


def lat_phi_grid( phi_min = -np.pi, phi_max = np.pi, lat_min = -np.pi / 2, lat_max = np.pi / 2, nlat = 180, nphi = 360 ):
    '''
    set up a 2d grid of latitude vs longitude
    with dimensions that will work well with imshow() later
    Note: the arrays are set up with dimensions [nlat, nphi] not the reverse, because of the way imshow works (it expects 2d array of [columns, rows])
    '''
    lat_1d = np.linspace( lat_min, lat_max, nlat)
    phi_1d = np.linspace( phi_min, phi_max, nphi)
    lat = ( (np.zeros([nphi, 1]) + 1.) @ lat_1d.reshape([1, nlat]) ).transpose()
    phi = (np.zeros([nlat, 1]) + 1.) @ phi_1d.reshape([1, nphi])
    return lat, phi
 

def get_all_cell_coordinates(f):
    # I believe these are volumetric (center of cells)
    cell_ids = f.read_variable('cellID')
    x = cell_ids * 0. ; y = cell_ids * 0. ; z = cell_ids * 0.
    for i, cell_id in enumerate(cell_ids):
        vec = f.get_cell_coordinates(cell_id)
        x[i] = vec[0]; y[i] = vec[1]; z[i] = vec[2]
    return x, y, z


def fg_grid(f, fg_b = None):
    # get the face positions of the cells
    # note that e.g. xmin, xmax give the leftmost and rightmost cell BOUNDARIES (the volumes are between these boundaries)
    dx = (f.read_parameter('xmax') - f.read_parameter('xmin')) / fg_b.shape[0]
    dy = (f.read_parameter('ymax') - f.read_parameter('ymin')) / fg_b.shape[1]
    dz = (f.read_parameter('zmax') - f.read_parameter('zmin')) / fg_b.shape[2]
    if (fg_b is None):
        fg_b = f.read_fsgrid_variable('fg_b')        # EGL: fg_b.shape = (1024, 736, 736, 3)
    x = np.linspace( f.read_parameter('xmin') + dx/2, f.read_parameter('xmax') - dx/2, fg_b.shape[0] )
    y = np.linspace( f.read_parameter('ymin') + dy/2, f.read_parameter('ymax') - dy/2, fg_b.shape[1] )
    z = np.linspace( f.read_parameter('zmin') + dz/2, f.read_parameter('zmax') - dz/2, fg_b.shape[2] )
    #x = np.linspace( f.read_parameter('xmin'), f.read_parameter('xmax'), fg_b.shape[0] )
    #y = np.linspace( f.read_parameter('ymin'), f.read_parameter('ymax'), fg_b.shape[1] )
    #z = np.linspace( f.read_parameter('zmin'), f.read_parameter('zmax'), fg_b.shape[2] )
    return x, y, z


def f_shue_parametrized(theta_polar, r_0, alpha):
    ''' Shue et al. (1997): A new functional form to study the solar wind control
    returns r as a function of theta,
    parameters: r_0, alpha ()
    '''
    return r_0 * (2 / (1 + np.cos(theta_polar)))**alpha

def f_shue_parameters(run):
    m_p  = 1.67262158e-27              # proton mass [kg]
    if (run == 'EGI'):
        B_z = -5                  # nT
        n_p = 1                   # cm^-3
        v_sw = 750                # km/sec
    elif (run == 'EGK'):
        B_z = -20                  # nT
        n_p = 1                   # cm^-3
        v_sw = 750                # km/sec
    elif (run == 'EGL'):
        B_z = -10
        n_p = 4
        v_sw = 750
    elif (run == 'EGM'):
        B_z = -5
        n_p = 1
        v_sw = 750
    elif (run == 'EGP'):
        B_z = -20    # B_x = -0.5??
        n_p = 7
        v_sw = 1000
    else:
        B_z = 0
        n_p = 1
        v_sw = 500
        print('VLASIATOR RUN NOT SPECIFIED!!!')     # error message
    n_p_SI = n_p * 1e6             # m^-3
    v_sw_SI = v_sw * 1e3           # m/sec
    rho = n_p_SI * m_p
    D_p = (rho / 2) * (v_sw_SI**2) * 1e9  # dynamic pressure, in nanoPAscals
    if B_z >= 0:
        r_0 = (11.4 + 0.013 * B_z) * D_p**(-1 / 6.6)              # Eq. 12, Shue 1997
    elif B_z < 0:
        r_0 = (11.4 + 0.14 * B_z) * D_p**(-1 / 6.6)
    alpha = ( 0.58 - (0.010 * B_z) ) * (1 + (0.010 * D_p))         # Eq. 13, Shue 1997
    return r_0, alpha


def f_shue(theta_polar, run = None):
    # Shue et al. (1997): 'A new functional form to study the solar wind control
    # of the magnetopause size and shape'
    # INPUTS 0 < theta_polar < 2pi
    # theta_polar is a numpy array
    # outputs magnetopause position r(theta) [in R_E]
    r_0, alpha = f_shue_parameters(run)
    r = f_shue_parametrized(theta_polar, r_0, alpha)
    return r, r_0, alpha


def f_shue_make_func(r_0, alpha):
    return lambda theta_polar: r_0 * (2 / (1 + np.cos(theta_polar)))**alpha



def interpolator_list_3d(f, vec_variable, bounds_error = True, fill_value = None):
    # vec_variable is a [nx,ny,nz,3] array
    # Create x, y, and z coordinates:
    # Read cellids in order to sort variables
    #cellids = vlsvReader.read_variable("CellID")
    xsize = vec_variable.shape[0]
    ysize = vec_variable.shape[1]
    zsize = vec_variable.shape[2]
    xmin = f.read_parameter('xmin')
    xmax = f.read_parameter('xmax')
    ymin = f.read_parameter('ymin')
    ymax = f.read_parameter('ymax')
    zmin = f.read_parameter('zmin')
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
        :returns:                  List of coordinates    
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
    tf_r_gt = (x_out**2 + y_out**2 + z_out**2)**0.5 > r_trace
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
        tf_r_gt_new = r_temp > r_trace
        mask = np.logical_and( np.logical_or( np.logical_and( ~tf_r_gt_new, tf_r_gt ), np.logical_and( tf_r_gt_new, ~tf_r_gt ) ), 
                              ~crossed_r_trace )      # test for which lines crossed r_trace this iteration, 
                                                      # requiring that r_trace has not been crossed before
        crossed_r_trace[mask] = True
        x_out[mask] = new_point[mask, 0]                # save where the field line crossed radius r_trace
        y_out[mask] = new_point[mask, 1]
        z_out[mask] = new_point[mask, 2]
        tf_r_gt = tf_r_gt_new    
        point = new_point
#        points.append( point )
    x_out[~crossed_r_trace] = point[~crossed_r_trace, 0]     # for remaining points that never cross the threshold, save their final location
    y_out[~crossed_r_trace] = point[~crossed_r_trace, 1]
    z_out[~crossed_r_trace] = point[~crossed_r_trace, 2]
    x_out = x_out.reshape(coord_list[0].shape)
    y_out = y_out.reshape(coord_list[0].shape)
    z_out = z_out.reshape(coord_list[0].shape)
    crossed_r_trace = crossed_r_trace.reshape(coord_list[0].shape)
    #    return points 
    return x_out, y_out, z_out, crossed_r_trace
    

def open_closed_boundary(bool_array):
    # Assumed 1D input Boolean array, works with slices of the output of test_field_open()
    # returns indices where of the array where the topology changes
    last_open_index_list = []
    last_closed_index_list = []
    for i in range(bool_array.size):
        if i > 0:
           if bool_array[i] == True and bool_array[i-1] == False:
               last_open_index_list.append(i)
               last_closed_index_list.append(i-1)
           if bool_array[i] == False and bool_array[i-1] == True:
               last_open_index_list.append(i-1)
               last_closed_index_list.append(i)
    return last_open_index_list, last_closed_index_list



def test_radial_field(f, coord_list_cart):
    # test whether the field lines are pointed away from the earth by tracing the line over a short distance (1 iteration)
    # returns a boolean array that can be used as a mask (True values means the field is quasi-radial)
    r0 = (coord_list_cart[0]**2 + coord_list_cart[1]**2 + coord_list_cart[2]**2)**0.5
    x, y, z, ind = trace_static_field( f, coord_list_cart, 1, 1, direction='+', fg_b = None, r_trace = 5 * R_EARTH )
    r = (x**2 + y**2 + z**2)**0.5
    return r > r0
    

def test_field_open(f, coord_list, fg_b = None, max_iterations = 15000, dx = R_EARTH / 50, trace_method = 'integrateB', coord_in='cartesian' ): 
    # return an array matching the elements of coord_list. True where the field lines are open (they leave the simulation), False where closed
    xlim = np.max([np.abs(f.read_parameter('xmin')), f.read_parameter('xmax')])
    ylim = np.max([np.abs(f.read_parameter('ymin')), f.read_parameter('ymax')])
    zlim = np.max([np.abs(f.read_parameter('zmin')), f.read_parameter('zmax')])
    rlim = (xlim**2 + ylim**2 + zlim**2)**0.5
    xpos, ypos, zpos, indpos = trace_static_field( f, coord_list, max_iterations, dx, direction='+', fg_b = fg_b, r_trace = rlim*1.1 )  #check the field line in both directions
    xneg, yneg, zneg, indneg = trace_static_field( f, coord_list, max_iterations, dx, direction='-', fg_b = fg_b, r_trace = rlim*1.1 )
    tf_nan = np.isnan(xpos * xneg * ypos * yneg * zpos * zneg)                     # trace_static_field returns nan when the trace leaves the simulation box 
    return tf_nan

 
def trace_coord(f, coord_list, fg_b = None, trace_method='dipole', direction = '+', coord_in='spherical', coord_out='cartesian', r_trace = 5 * R_EARTH, max_iterations = 1000, dx = R_EARTH / 50):
    '''
    trace_coord() takes input coordinates given by field_list and maps them to a location in the magnetosphere by following trajectories along field lines.
    
    f: a VlsvReader object, from reading a .vlsv file
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
    '''
    elif trace_method == 'particlePusher':
       x_out, y_out, z_out = np.zeros(len(x)), np.zeros(len(x)), np.zeros(len(x))   # dummy
       # *** fill in this code, generate x_out, y_out, z_out
    '''
    if coord_out == 'spherical':
        r_out, theta_out, phi_out = cartesian_to_spherical(x_out, y_out, z_out)
        coord_list_out = [r_out, theta_out, phi_out]
    elif coord_out == 'cartesian':
        coord_list_out = [x_out, y_out, z_out]
    return coord_list_out      

 
def get_variable(f, coord_list_cart, var_list = ['proton/vg_precipitationdifferentialflux'], operator = "pass", interpolation_order = 1):
    ''' get_variable() takes input coordinates specified by coord_list ( e.g., the output of trace_coord) and computes the value of a grid variable at that location
 
    	INPUTS:   
    	:param f:			a VlsvReader object, from reading a .vlsv file
    	:param coord_list_cart:		a 3-element list of 1D numpy arrays of equal length N (or scalars?),
					specifying the coordinates (GSE xyz) to generate variables in var_list
					ex. [ [1.2e7 2.3e7 3.1e7] , [0 0 0] , [1e7 2e7 3e7] ]
					Note: nans are NOT accepted in elements of coord_list_cart, need to handle this outside of the function
	:kword var_list:		a list of the variables to be returned at the specified N points
 	:kword operator: 		see vlsv_intpol_points() in pyCalculations module
   	:kword interpolation_order:	self-explanatory, but see vlsv_intpol_points() in pyCalculations module. NOTE: testing interpolation_order=2 led to a lot of zeros...
	:returns:			a tuple (coord_out, cellids, params)
    		:coord_out: 		(N, 3) numpy array of coordinates (this is just identical to contents of coord_list_cart?)
		:cellids:		cellids used for interpolation
		:params:		a dictionary with X entries (each is a N-element numpy array), where X is ithe number of variables in var_list_output. Note len(var_list_output) >= len(var_list)   
    '''
    points = (np.array([coord_list_cart[0], coord_list_cart[1], coord_list_cart[2]] )).transpose()
    coord_out, cellids, params, hstr = vlsv_intpol_points(f, points, var_list, operator=operator, interpolation_order=interpolation_order)
    var_list_output = hstr.split()[4:]
    dct = {}
    for i, varname in enumerate(var_list_output):
        dct[varname] = params[:,i]
    return coord_out, cellids, dct


def volume_shell(data, x, y, z, threshold, delta, step = 1): #STEP = 1 OR STEP = 0.5
    ''' Consider data that are within delta of the threshold
        Then loop through rcyl (cylindrical polar coordinate), and 
        for each iteration select the point with the greatest x.
        (note solar wind comes from +x direction). Combine all
        these data to define the shell
    '''
    x_min = np.min(np.abs(x))
    x_index = np.where(x>=x_min)[0]
    new_z = z[x_index]
    new_y = y[x_index]
    new_x = x[x_index]
    new_data = data[x_index]

    z_min = np.min(np.abs(new_z)) # ~ 0.08 R_E (500 km)
    z_index = np.where(np.abs(new_z)<=0.5)[0] # where |z|<=0.5 R_E
    new_z = new_z[z_index]
    new_y = new_y[z_index]
    new_x = new_x[z_index]
    new_data = new_data[z_index]

    rcyl = np.sign(new_y) * (new_y**2 + new_z**2)**0.5 #MAGNETOPAUSE RADIAL DISTANCE DEPENDS ON Y & Z
    rcyls = np.arange(-5, 5, step)
    #rcyls = np.arange(np.nanmin(rcyl), np.nanmax(rcyl), step) #list of sqrt(y^2+z^2) data points
    x_shell = []; rcyl_shell = []
    for rcyl_i in rcyls:
        ind1, = np.nonzero((new_data >= threshold) & (rcyl>=rcyl_i) & (rcyl<(rcyl_i+step)))
        #ind1, = np.nonzero((np.abs(new_data - threshold) < delta) & (rcyl>=rcyl_i) & (rcyl<(rcyl_i+step)))
        
        #ind1, = np.nonzero(((new_data - threshold) > 0) & (rcyl>=rcyl_i) & (rcyl<(rcyl_i+step)))   # NOTE: this deprecates delta
        if ind1.size > 0:
            x_shell.append(np.nanmin(new_x[ind1]))
            rcyl_shell.append(rcyl_i + step/2)
    return np.array(x_shell), np.array(rcyl_shell)
    #return np.abs(new_data - threshold) < delta


def bow_shock(f, fg_b = None, threshold = None, delta = 0.05, method = 'n'):
    ''' define bow shock as contour where some plasma parameter is equal to some threshold
        TODO: IMPLEMENT A BETTER MODEL (for now, using Shue)
    '''
    vg_x, vg_y, vg_z = get_all_cell_coordinates(f)
    x = vg_x / R_EARTH; y = vg_y / R_EARTH; z = vg_z / R_EARTH
    #vg_v = f.read_variable('vg_v')
    #vg_vtot = (vg_v[:,0]**2 + vg_v[:,1]**2 + vg_v[:,2]**2)**0.5
    if method == 'n':
        ''' following Battarbee (2020), use density contour to define large-scale bow shape
        '''
        if threshold == None:
            threshold = 2
        xmax_list = [ np.array([f.read_parameter('xmax')-1, 0, 0])]     # position @ xmax boundary (the -1 is to make sure it falls within a cell)
        coord_out, cellids, params, hstr = vlsv_intpol_points(f, xmax_list, ['proton/vg_rho'], operator='pass', interpolation_order=1)
        n_sw = params[0].flatten()
        check = f.read_variable('proton/vg_rho') / n_sw
    elif method == 'mach_ma':
        # Alfvenic mach number
        if threshold == None:
            threshold = 1
        check = f.read_variable('vg_ma')     # v5 data reducer for Mach number, see reduction.py
    elif method == 'mach_ms':
        # magnetosonic mach number
        if threshold == None:
            threshold = 1
        check = f.read_variable('vg_mms')     # v5 data reducer for Mach number, see reduction.py
    elif method == 'b':
        ''' check for where a significant change in B is observed, relative to xmax plane
            For now, just look for changes in the magnitude of B. But dot product of unit vectors is probably better
           Note: will need to re-implement t
        '''
        if threshold == None:
            threshold = 1
        #if (fg_b is None):
        #    fg_b = f.read_variable('fg_b')
        #fg_x, fg_y, fg_z = fg_grid(f, fg_b = fg_b)
        B3 = f.read_variable('vg_b_vol')
        B = (B3[:,0]**2 + B3[:,1]**2 + B3[:,2]**2)**0.5
        #B = (fg_b[:,:,:,0]**2 + fg_b[:,:,:,1]**2 + fg_b[:,:,:,2]**2)**0.5
        #b_hat = fg_b
        b_hat = B3
        #b_hat[:,:,:,0] = b_hat[:,:,:,0]/B; b_hat[:,:,:,1] = b_hat[:,:,:,1]/B; b_hat[:,:,:,2] = b_hat[:,:,:,2]/B 
        b_hat[:,0] = b_hat[:,0]/B; b_hat[:,1] = b_hat[:,1]/B; b_hat[:,2] = b_hat[:,2]/B 
        #b_hat_ref = b_hat[-1, int(b_hat.shape[1]/2), int(b_hat.shape[2]/2), :]           #x = xmax, y = 0, z = 0
        b_hat_ref = np.array([0, 0, -1])          #x = xmax, y = 0, z = 0
        #dot = b_hat[:,:,:,0]*b_hat_ref[0] + b_hat[:,:,:,1]*b_hat_ref[1] + b_hat[:,:,:,2]*b_hat_ref[2]
        check = b_hat[:,0]*b_hat_ref[0] + b_hat[:,1]*b_hat_ref[1] + b_hat[:,2]*b_hat_ref[2]    # projection
    x_shell, rcyl_shell = volume_shell(check, x, y, z, threshold, delta)
    return x_shell, rcyl_shell




def magnetopause(f, threshold = 1, delta = 0.05):
    ''' define magnetospause as where modified plasma beta equals 1
        (Brenner et al, 2021)
    '''
    vg_x, vg_y, vg_z = get_all_cell_coordinates(f)
    x = vg_x / R_EARTH; y = vg_y / R_EARTH; z = vg_z / R_EARTH
    #vg_v = f.read_variable('vg_v')
    #vg_vtot = (vg_v[:,0]**2 + vg_v[:,1]**2 + vg_v[:,2]**2)**0.5
    P_th = f.read_variable('proton/vg_pressure') 
    P_dyn = f.read_variable('proton/vg_pdyn')
    B3 = f.read_variable('vg_b_vol')
    B = (B3[:,0]**2 + B3[:,1]**2 + B3[:,2]**2)**0.5
    beta_prime = (P_th + P_dyn) / ( B**2 / (2 * mu_0) )
    x_shell, rcyl_shell = volume_shell(beta_prime, x, y, z, threshold, delta)
    return x_shell, rcyl_shell
    #return x[mask], y[mask], z[mask]



def fit_bow_shock(f, run = '', root_dir = '', fileIndex = '', threshold = 2, delta = 0.05, plot = True, method = 'n'):
    #vg_x, vg_y, vg_z = bow_shock(f, threshold = threshold, delta = delta)
    #vg_rcyl = (vg_y**2 + vg_z**2)**0.5
    vg_x, vg_rcyl = bow_shock(f, threshold = threshold, delta = delta, method = method)
    r_bs = (vg_rcyl**2 + vg_x**2)**0.5
    theta_bs = np.arccos( vg_x / r_bs )
    r_0, alpha = 10, 0.6
    params, params_cov = scipy.optimize.curve_fit(f_shue_parametrized, theta_bs, r_bs, p0=[r_0, alpha])
    f_shue_fit = f_shue_make_func(params[0], params[1])
    #PLOT 1: x vs r_cyl of bow shock points
    if plot:
        save_dir = '{}{}/magnetopause_pos/{}/'.format(root_dir, run.upper(), str(fileIndex).zfill(7))         #CUSTOM
        rplot = f_shue_fit(theta_bs)
        title = r'{} Bow shock position, time={}'.format(run.upper(), fileIndex)
        plt.scatter(vg_x, vg_rcyl, color = 'blue')
        plt.xlabel(r'x [$r_E$]')
        plt.ylabel(r'$\sqrt{y^2 + z^2}$ [$r_E$]')
        plt.title(title)
        plt.plot(rplot * np.cos(theta_bs), rplot * np.sin(theta_bs))
        filename = '{}bow_shock_pos_data_{}_{}.png'.format(save_dir, method, str(fileIndex).zfill(5))
        plt.savefig(filename)
        plt.close()
    return f_shue_fit, theta_bs, r_bs


def fit_magnetopause(f, run = '', root_dir = '', fileIndex = '', threshold = 1, delta = 0.05, plot = True):
    #vg_x, vg_y, vg_z = bow_shock(f, threshold = threshold, delta = delta)
    #vg_rcyl = (vg_y**2 + vg_z**2)**0.5
    vg_x, vg_rcyl = magnetopause(f, threshold = threshold, delta = delta)
    r_mp = (vg_rcyl**2 + vg_x**2)**0.5
    theta_mp = np.arccos( vg_x / r_mp )
    r_0, alpha = 10, 0.6
    params, params_cov = scipy.optimize.curve_fit(f_shue_parametrized, theta_mp, r_mp, p0=[r_0, alpha])
    f_shue_fit = f_shue_make_func(params[0], params[1])
    #PLOT 1: x vs r_cyl of magnetopause points
    if True:
        save_dir = root_dir       #CUSTOM
        #---------------------------------------------------------------------------------
        #run = 'EGL'
        #dim = '3D'
        #bulk = 'bulk1.egl'

        # Defining source and output file locations
        #bulkLocation = '/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/'.format(dim, run)
        outputLocation = root_dir

        # Defining source and output file locations
        run = 'EGP'   # EGP
        outputLocation = '/wrk-vakka/users/horakons/carrington/plots/{}_validation_paper/'.format(run)
        if run == 'EGL':
            j=857                       #EGL --- pressure pulse arrival time
            before_index = 621 #     EGL: 621
            after_index = 1760 #      EGL: 1760
        elif run == 'EGP':
            j=506                       #EGP --- first time with precipitation data   (269-506 total)
            before_index = 352 #     EGL: 621
            after_index = 506 #      EGL: 1760

        # Time step number
        j = int(fileIndex)
        # Bulk file name
        #bulkname1 = '{}.{}.vlsv'.format(bulk,str(j).zfill(7))
        bulkname = get_vlsvfile_fullpath(run, j)
       
        #---------------------------------------------------------------------------------
       
        ax = plt.gca()

        #Shue (1998):
        if j <= 857: #before shock
            D = 1e6 * 1.67262192e-27 * (7.5e5)**2 * 1e9 #nPa
            B = -5 #nT
        else: #after shock
            D = 4e6 * 1.67262192e-27 * (7.5e5)**2 * 1e9 #nPa
            B = -10 #nT

        alpha = (0.58 - 0.007 * B) * (1 + 0.024 * np.log(D))
        r_0 = (10.22 + 1.29 * np.tanh(0.184 * (B + 8.14))) * D ** (-1/6.6)
        theta = np.linspace(0, 2*np.pi, 1000)
        r = r_0 * (2/(1 + np.cos(theta))) ** alpha
        
        m_x = r*np.cos(theta)
        m_R = r*np.sin(theta)
        ax.plot(m_x, m_R, color='cyan', linewidth=2.5, zorder=1)

        #beta_star colormap:
        streamlines = 'proton/vg_v'
        sl_color = 'orange'
        sl_density = 2
        normal = 'z' #y=0 or z=0 plane
        scale = 1.5
       
            #curve fit: 
        model = np.poly1d(np.polyfit(vg_rcyl, vg_x, 2))
        polyline = np.linspace(-15, 15, 1000) #1000 data points
        #subsolar point:
        mp_nose_index = np.where(polyline == np.min(np.abs(polyline)))
        mp_nose_x = model(polyline)[mp_nose_index][0]
        mp_nose_y = polyline[mp_nose_index]

        if plot:
            #pt.plot.plot_colormap3dslice(filename=bulkLocation+bulkname1,var='proton/vg_beta_star', boxre=[-15, 15, -15, 15], normal =normal, run=run,
            pt.plot.plot_colormap3dslice(filename=bulkname,var='proton/vg_beta_star', boxre=[-15, 15, -15, 15], normal =normal, run=run,
                                         colormap='seismic',vmin=10**(-1)*0.5,vmax=10*0.5,step=j,outputdir=outputLocation,
                                         outputfile='beta_star_colormap_xy_{}_{}.pdf'.format(run, str(j).zfill(5)), 
                                         Earth=1, streamlines=streamlines, streamlinedensity=sl_density, streamlinecolor = sl_color, cutpointre=0, axes = ax, scale=scale, useimshow=True)
       
            mp_col = 'lime'  # magenta
            ax.plot(model(polyline), polyline, color=mp_col, linewidth=2.5, zorder=3)

            #beta_star circles:
            ax.scatter(vg_x, vg_rcyl, color='w', zorder=2) #plot in xy-plane
        
            ax.plot(mp_nose_x, mp_nose_y, marker='o', mfc=mp_col, mec=mp_col, zorder=4)
        
            #----------------------------------------------------------------------------
            filename = '{}beta_star_colormap_xy_{}.pdf'.format(save_dir, str(fileIndex).zfill(5))
            plt.savefig(filename)
            plt.close()


            ax = plt.gca()

            normal = 'y'
            streamlines = 'vg_b_vol'
            sl_color = 'yellow'
            ax.plot(m_x, m_R, color='cyan', linewidth=2.5, zorder=1)

            #pt.plot.plot_colormap3dslice(filename=bulkLocation+bulkname1,var='proton/vg_beta_star', boxre=[-15, 15, -15, 15], normal =normal, run=run,
            pt.plot.plot_colormap3dslice(filename=bulkname,var='proton/vg_beta_star', boxre=[-15, 15, -15, 15], normal =normal, run=run,
                                         colormap='seismic',vmin=10**(-1)*0.5,vmax=10*0.5,step=j,outputdir=outputLocation,
                                         outputfile='beta_star_colormap_xz_{}_{}.pdf'.format(run, str(j).zfill(5)), 
                                         Earth=1, streamlines=streamlines, streamlinedensity=sl_density, streamlinecolor = sl_color, cutpointre=0, axes = ax, scale=scale, useimshow=True)
            #beta_star circles:
            #ax.scatter(vg_x, vg_rcyl, color='w', zorder=2) #plot in xy-plane

        #curve fit: 
        model = np.poly1d(np.polyfit(vg_rcyl, vg_x, 2))
        polyline = np.linspace(-15, 15, 1000) #1000 data points
        #subsolar point:
        mp_nose_index = np.where(polyline == np.min(np.abs(polyline)))
        mp_nose_x = model(polyline)[mp_nose_index][0]
        mp_nose_y = polyline[mp_nose_index]

        if plot:
            ax.plot(model(polyline), polyline, color=mp_col, linewidth=2.5, zorder=3)
            ax.plot(mp_nose_x, mp_nose_y, marker='o', mfc=mp_col, mec=mp_col, zorder=4)
            #----------------------------------------------------------------------------
            filename = '{}beta_star_colormap_xz_{}_{}.pdf'.format(save_dir, run, str(fileIndex).zfill(5))
            plt.savefig(filename)
            plt.close()

        print(j, mp_nose_x)
        return mp_nose_x


def label_region_func(f):
    ''' returns a function that can take cartesian coordinates as an argument
        and outputs an array of labels
    '''
    f_mp, theta_mp, r_mp = fit_magnetopause(f, plot = False)
    f_bs, theta_bs, r_bs = fit_bow_shock(f, plot = False)
    vg_x, vg_y, vg_z = get_all_cell_coordinates(f)
    vg_r_min = np.nanmin( (vg_x**2 + vg_y**2 + vg_z**2)**0.5 )
    def output_func(x, y, z):
        ''' x, y, z are numpy arrays of equal shape
            output a character array that specifies the region:
            'sw': Solar Wind
            'sh': magnetoSHeath
            'mg': MaGnetosphere --- anything behind the magnetopause
            'ob': Out-of-Bounds --- anything outside of Vlasiator simulation
            Note that these regions are non-overlapping
        '''
        if type(x) == int or type(x) == float:
            return output_func( np.array([x]), np.array([y]), np.array([z]) )[0]
        r = (x**2 + y**2 + z**2)**0.5
        theta = np.arccos( x / r )
        labels = np.char.array(x, unicode = False, itemsize = 2)
        labels[:] = 'mg'
        ind_sw = r > f_bs(theta)
        labels[ind_sw] = 'sw'
        ind_sh = (r <= f_bs(theta)) & (r>f_mp(theta))
        labels[ind_sh] = 'sh'
        ind_ob = (x > f.read_parameter('xmax')) | (x < f.read_parameter('xmin')) | \
                 (y > f.read_parameter('ymax')) | (y < f.read_parameter('ymin')) | \
                 (z > f.read_parameter('zmax')) | (z < f.read_parameter('zmin')) | \
                 (r < np.nanmin(vg_r_min))
        labels[ind_ob] = 'ob'         
        return labels
    return output_func



def label_region(f, x, y, z):
    ''' x, y, z are numpy arrays of equal shape
        output a character array that specifies the region:
        'sw': Solar Wind
        'sh': magnetoSHeath
        'mg': MaGnetosphere --- anything behind the magnetopause
        'ob': Out-of-Bounds --- anything outside of Vlasiator simulation
        Note that these regions are non-overlapping
    '''
    return label_region_func(f)(x, y, z)



class ParamObj:
    ''' The ParamObj class plots figures, returns data, and saves that data for a given data frame
        Used for analyzing the >=5 data signatures associated with solar storms:
        dB/dt, FAC, proton flux, polar cap boundary position, magnetopause position...
        Each of these parameters will have its own subclass of ParamObj, so that code doesn't have to be repeated

    :param f:			a VlsvReader object, from reading a .vlsv file
    :param run:			the run name, e.g. 'EGL'
    :param coord_list:		a 3-element list of numpy arrays of equal length (or scalars?),
              			specifying the GSE [SI] position, spherical or cartesian as designated by coord_in
				ex. [ [1.2e7 2.3e7 3.1e7] , [0 0 0] , [1e7 2e7 3e7] ]
    '''
    def __init__(self, run, f, coord_list, coord_list_traced = None, fileIndex = None, file_suffix = '', plot_data=False, make_plot_data = False, coord_in = 'spherical'):
        self.root_dir = ROOT_DIR
        self.run = run
        self.f = f
        if coord_in == 'cartesian':
            self.coord_list = coord_list
            self.coord_list_sphere = list(cartesian_to_spherical(*self.coord_list))
        elif coord_in == 'spherical':
            self.coord_list_sphere = coord_list
            self.coord_list = list(spherical_to_cartesian(*self.coord_list_sphere))
        if coord_list_traced is None:
            self.coord_list_traced = self.coord_list
            self.file_suffix = '' + file_suffix
        else:
            self.file_suffix = '_traced' + file_suffix
            if coord_in == 'cartesian':
                self.coord_list_traced = coord_list_traced
            elif coord_in == 'spherical':
                self.coord_list_traced =list(spherical_to_cartesian( *coord_list_traced ))
        if fileIndex is None:        #klug
           self.fileIndex = f.read_parameter('fileIndex')
        else:
           self.fileIndex = fileIndex
        self.shape = coord_list[0].shape
        self.size = coord_list[0].size
        self.data_list = []          # data associated with particular coordinates
        self.data_label_list = []    # ""
        self.fig_path_list = []      # list of paths where figures will be saved
        self.fig_list = []                   # list of figures (matplotlib plt objects)
        self.coord_mask = ~np.isnan( self.coord_list[0] * self.coord_list[1] * self.coord_list[2] *
                                     self.coord_list_traced[0] * self.coord_list_traced[1] * self.coord_list_traced[2] )          # test for where coordinates are finite
        if make_plot_data:
            self.make_plot_data(plot_data=plot_data)
    def merge(self, pobjs):
        if type(pobjs) == list:
            for p in pobjs:
                self.merge(p)
        else:
            self.fig_list.extend(pobjs.fig_list)
            self.fig_path_list.extend(pobjs.fig_path_list)
            self.data_list.extend(pobjs.data_list)
            self.data_label_list.extend(pobjs.data_label_list)
    def get_variable(self, var_list = None):
        coord_list_gv = [ self.coord_list_traced[0][self.coord_mask],
                          self.coord_list_traced[1][self.coord_mask],
                          self.coord_list_traced[2][self.coord_mask] ]             # coordinates to be input into get_variable()
        coord_out, cellids, dct = get_variable(self.f, coord_list_gv, var_list = var_list)
        return coord_out, cellids, dct
#    def mkdir_path(self, path):
#        filedir_list = path.split('/')
#        filedir = path[:-len(filedir_list[-1])]
#        if not(os.path.exists(filedir)):
#            os.system('mkdir -p {}'.format(filedir))                  # need a function to just compute the save dir from the path using split('/') and indexing
    def make_lonlat_plot(self, phi, theta, plot_variable_2d, xlabel='longitude [deg.]', ylabel='latitude [deg.]', title = '', cbar = True, cbar_label='[SI]', extent = [-180, 180, -90, 90], 
                         vmin = None, vmax = None, cmap = plt.cm.get_cmap('plasma'), path = './test.png'):
        # PLOT 1: azimuthal projection of fluxes
        fig, ax = plt.subplots()
        #im = ax.imshow(plot_variable_2d, extent = extent, interpolation = None,
        #               cmap = cmap, vmin = vmin, vmax = vmax)
        im = ax.pcolormesh(phi*180/np.pi, theta*180/np.pi, plot_variable_2d,
                           cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if cbar:
            divider = make_axes_locatable(ax)
            cax = divider.append_axes('right', size='5%', pad=0.05)
            mycbar = fig.colorbar(im, cax=cax, orientation='vertical')
            #mycbar = ax.figure.colorbar(im)
            mycbar.set_label(cbar_label)
        #plt_temp = plt
        # self.fig_list.append(plt_temp)
        #self.fig_path_list.append(path)
        plt.tight_layout()
        mkdir_path(path)
        plt.savefig(path) 
        plt.close()
        # add plt, filename to fig_list and fig_path_list, respectively (or implement as a dictionary instead?)
    def make_polar_plot(self, phi, theta, plot_variable_2d, cbar = True, cbar_label='[SI]', title = '',
                         vmin = None, vmax = None, cmap = plt.cm.get_cmap('plasma'), path = './test.png'):
        fig = plt.figure()
        gs = gridspec.GridSpec(1, 2, width_ratios=[10,1])
        #gs = gridspec.GridSpec(1, 1)
        ax1 = plt.subplot(gs[0], projection="polar", aspect=1.)
        phi_plot = phi + (np.pi / 2)          # rotate by ninety-degrees so noon is at the top
        theta_plot = theta*180/np.pi
        im = ax1.pcolormesh(phi_plot, theta_plot, plot_variable_2d,
                            cmap=cmap, vmin=vmin, vmax=vmax, shading='auto')            # north pole
        ax1.set_ylim([0,35])
        ax1.set_yticks([10, 20, 30])  # Less radial ticks
        ax1.set_yticklabels([r'$80^\circ$', r'$70^\circ$', r'$60^\circ$'])
        ax1.set_rlabel_position(-10)  # Move radial labels away from plotted line
        ax1.set_xticks(list(np.arange(0, 2*np.pi, np.pi / 4)))
        #ax1.set_xticklabels(list(np.arange(18, 24, 3)) + list(np.arange(0, 18, 3)))
        ax1.set_xticklabels(list(np.arange(6, 24, 3)) + list(np.arange(0, 6, 3)))
        ax1.grid(True)
        if cbar:
           ax2 = plt.subplot(gs[1])
           #cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
           #plt.colorbar(im, cax=cax) # Similar to fig.colorbar(im, cax = cax)
           #divider = make_axes_locatable(ax1)
           #cax = divider.append_axes('right', size='5%', pad=0.05)
           #cbar = fig.colorbar(im, cax=cax, orientation='vertical')
           cbar = plt.colorbar(im, cax=ax2, aspect = 40)
           #cbar = mpl.colorbar.ColorbarBase(ax2, cmap = cmap)
           cbar.set_label(cbar_label)
        ax1.set_title(title)
        #plt_temp = plt
        #self.fig_list.append(plt_temp)
        #self.fig_path_list.append(path)
        plt.tight_layout()
        mkdir_path(path)
        plt.savefig(path) 
        plt.close()
    def write_data(self, filename='data.csv', coord_in='spherical', append=False):
        dct = {}
        if coord_in == 'cartesian':
            dct['x_GSE_m'] = list(self.coord_list[0].flatten())
            dct['y_GSE_m'] = list(self.coord_list[1].flatten())
            dct['z_GSE_m'] = list(self.coord_list[2].flatten())
            dct['x_GSE_traced_m'] = list(self.coord_list_traced[0].flatten())
            dct['y_GSE_traced_m'] = list(self.coord_list_traced[1].flatten())
            dct['z_GSE_traced_m'] = list(self.coord_list_traced[2].flatten())
        elif coord_in == 'spherical':
            dct['r_GSE_m'] = list(self.coord_list_sphere[0].flatten())
            #dct['theta_GSE_radians'] = list(self.coord_list_sphere[1].flatten())
            #dct['phi_GSE_radians'] = list(self.coord_list_sphere[2].flatten())
            dct['theta_GSE_deg'] = list(self.coord_list_sphere[1].flatten() * 180 / np.pi)
            dct['phi_GSE_deg'] = list(self.coord_list_sphere[2].flatten() * 180 / np.pi)
            r_traced, theta_traced, phi_traced = cartesian_to_spherical(self.coord_list_traced[0], self.coord_list_traced[1], self.coord_list_traced[2]) 
            dct['r_GSE_traced_m'] = list(r_traced.flatten())
            #dct['theta_GSE_traced_radians'] = list(theta_traced.flatten())
            #dct['phi_GSE_traced_radians'] = list(phi_traced.flatten())
            dct['theta_GSE_traced_deg'] = list(theta_traced.flatten() * 180 / np.pi)
            dct['phi_GSE_traced_deg'] = list(phi_traced.flatten() * 180 / np.pi)
        dct['t_sec'] = np.zeros(self.coord_list[0].size)  + self.f.read_parameter('time')
        for label, data in zip(self.data_label_list, self.data_list):
            dct[label] = list(data.flatten())
        df = pandas.DataFrame.from_dict(dct)
        if int(ARGS.nproc)  > 1:
            #when multi-threading,make many files that are all concatenated together 
            filename_csv = filename[0:-4] + '_' + str(self.fileIndex).zfill(4) + '.csv'
        else:
            filename_csv = filename
        mkdir_path(filename_csv)
        df_write = df*1.0  #converts True to 1.0, False to 0.0
        if append:
            df_write.to_csv(filename_csv, mode='a', header=not os.path.exists(filename), index=False)
        else:
            df_write.to_csv(filename_csv, index=False)
    def savefigs(self):
        # save the data used to make the figures, in a format that can be read by FMI
        for fig, path in zip(self.fig_list, self.fig_path_list):
            mkdir_path(path)
            fig.savefig(path)


class pFluxObj(ParamObj):
    def __init__(self, run, f, coord_list, **kwargs):
        ParamObj.__init__(self, run, f, coord_list, **kwargs)
    def make_plot_data(self, plot_data=False):
        r, theta, phi = self.coord_list_sphere
        coord_out, cellids, dct = self.get_variable(var_list = ['proton/vg_precipitationdifferentialflux'])     # store variables to plot in dictionary (dct)
        for i, key in enumerate(dct):
            proton_energy = self.f.read_parameter('proton_PrecipitationCentreEnergy{}'.format(i))
            title = '{} {} eV proton DEF [SI], time={}'.format(self.run.upper(), int(proton_energy), self.fileIndex)
            plot_variable_2d = np.zeros(self.size) * np.nan
            plot_variable_2d[self.coord_mask.reshape(self.size)] = dct[key]
            plot_variable_2d = plot_variable_2d.reshape(self.shape)
            self.data_list.append(plot_variable_2d)
            self.data_label_list.append( 'proton_DEF_{}_eV'.format(str(int(proton_energy)).zfill(5)) )
            # PLOT 1: azimuthal projection of fluxes
            delta_lat = 40                   # range of latitudes (in degrees) to plot
            nlat = self.shape[0]
            delta_lat_ind = int(nlat * (delta_lat / 180))
            save_dir = '{}{}/proton_DEF/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))
            filename = '{}azim_frame_{}_{}_ev{}.png'.format(save_dir, str(self.fileIndex).zfill(5), str(int(proton_energy)).zfill(5), self.file_suffix)
            #self.make_polar_plot(phi[(nlat-delta_lat_ind-1):,:], theta[(nlat-delta_lat_ind-1):,:], plot_variable_2d[(nlat-delta_lat_ind-1):,:], 
            #                     cbar_label=r'Diff. energy flux [cm$^{-2} s^{-1} sr^{-1} eV^{-1}$]', title = title,
            #                     vmin = 0, vmax = np.nanmax(plot_variable_2d), cmap = plt.cm.get_cmap('plasma'), path = filename)
            if plot_data:
                self.make_polar_plot(phi[(nlat-delta_lat_ind-1):,:], theta[(nlat-delta_lat_ind-1):,:], np.log10(plot_variable_2d[(nlat-delta_lat_ind-1):,:]),
                                     cbar_label=r'log$_{10}$(Diff. energy flux) [cm$^{-2} s^{-1} sr^{-1} eV^{-1}$]', title = title,
                                     #vmin = np.nanmin(np.log10(plot_variable_2d[(nlat-delta_lat_ind-1):,:] )), 
                                     #vmax = np.nanmax(np.log10(plot_variable_2d[(nlat-delta_lat_ind-1):,:] )), 
                                     vmin = np.log10(100000*np.exp(-proton_energy / 5000)) - 10, 
                                     vmax = np.log10(100000*np.exp(-proton_energy / 5000)), 
                                     cmap = plt.cm.get_cmap('plasma'), path = filename)
                # PLOT 2: square long-lat projection of fluxes
                save_dir = '{}{}/proton_DEF/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))
                filename = '{}rect_frame_{}_{}_ev{}.png'.format(save_dir, str(self.fileIndex).zfill(5), str(int(proton_energy)).zfill(5), self.file_suffix)
                self.make_lonlat_plot(phi, theta, plot_variable_2d, xlabel='longitude [deg.]', ylabel='latitude [deg.]', title = title, cbar_label='', extent = [-180, 180, -90, 90],
                                      vmin = 0, vmax = 10000*np.exp(-proton_energy / 5000), cmap = plt.cm.get_cmap('plasma'), path = filename)


class dBdtObj(ParamObj):
    def __init__(self, run, f, coord_list, **kwargs):
        self.dct = None
        ParamObj.__init__(self, run, f, coord_list, **kwargs)
    def make_plot_data(self, plot_Data=False):
        #f_2 = pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/{}/bulk/bulk1.{}.{}.vlsv".format(self.run.upper(), self.run.lower(), str(self.fileIndex-1).zfill(7) ))
        f_2 = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(self.run, self.fileIndex - 1))
        dBO_2 = dBdtObj(self.run, f_2, self.coord_list_sphere, make_plot_data = False )
        coord_out, cellids, dct = self.get_variable(var_list = ['vg_b_vol'])     # store variables to plot in dictionary (dct)
        coord_out_2, cellids_2, dct_2 = dBO_2.get_variable(var_list = ['vg_b_vol'])     # store variables to plot in dictionary (dct)
        self.dct = dct
        r, theta, phi = self.coord_list_sphere
        subscripts = ['x', 'y', 'z']
        for i, key in enumerate(dct):
            # PLOT(s): dB_x/dt, dB_y/dt, dB_z/dt
            plot_variable_2d = np.zeros(self.size) * np.nan
            plot_variable_2d[self.coord_mask.reshape(self.size)] = dct[key]
            plot_variable_2d = plot_variable_2d.reshape(self.shape)
            plot_variable_2d_2 = np.zeros(dBO_2.size) * np.nan
            plot_variable_2d_2[dBO_2.coord_mask.reshape(dBO_2.size)] = dct_2[key]
            plot_variable_2d_2 = plot_variable_2d_2.reshape(dBO_2.shape)
            plot_variable_2d = (plot_variable_2d - plot_variable_2d_2) / self.f.read_parameter('dt')             # take the time derivative!
            title = r'{}, $dB_{}/dt$, time={}'.format(self.run.upper(), subscripts[i], self.fileIndex)
            self.data_list.append(plot_variable_2d)
            self.data_label_list.append( '$dB{}_dt_SI$'.format(subscripts[i]) )
            if plot_data:
                save_dir = '{}{}/dBdt/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))         #CUSTOM
                filename = '{}rect_frame_{}_dB{}dt{}.png'.format(save_dir, str(self.fileIndex).zfill(5), subscripts[i], self.file_suffix)
                self.make_lonlat_plot(phi, theta, plot_variable_2d, xlabel='longitude [deg.]', ylabel='latitude [deg.]', title = title, cbar_label=r'dB$_{}$/dt [T/sec]'.format(subscripts[i]),
                                      extent = [-180, 180, -90, 90], vmin = np.nanmin(plot_variable_2d), vmax = np.nanmax(plot_variable_2d), cmap = plt.cm.get_cmap('plasma'), path = filename) 




class FieldOpenObj(ParamObj):
    def __init__(self, run, f, coord_list, fg_b = None, max_iterations = 15000, dx = R_EARTH / 50, **kwargs):
        if (fg_b is None):
            fg_b = f.read_variable('fg_b')     
        self.fg_b = fg_b
        self.max_iterations = max_iterations
        self.dx = dx
        ParamObj.__init__(self, run, f, coord_list, **kwargs) 
    def make_plot_data(self, plot_data=False):
        plot_variable_2d = test_field_open(self.f, self.coord_list, fg_b = self.fg_b, max_iterations = self.max_iterations, dx = self.dx, trace_method = 'integrateB', coord_in = 'cartesian' )
        self.data_list.append(plot_variable_2d)
        self.data_label_list.append('open_vs_closed')
        r, theta, phi = self.coord_list_sphere
        save_dir = '{}{}/field_open/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))         #CUSTOM
        if plot_data:
            #PLOT 0: 
            delta_lat = 40                   # range of latitudes (in degrees) to plot
            nlat = self.shape[0]
            delta_lat_ind = int(nlat * (delta_lat / 180))
            filename = '{}azim_frame_{}_fieldopen{}.png'.format(save_dir, str(self.fileIndex).zfill(5), self.file_suffix)
            title = r'{} Field line topology, time={}'.format(self.run.upper(), self.fileIndex)
            self.make_polar_plot(phi[(nlat-delta_lat_ind-1):,:], theta[(nlat-delta_lat_ind-1):,:], plot_variable_2d[(nlat-delta_lat_ind-1):,:],
                                 #cbar_label='Open [1] vs. closed regions [0]',
                                 cbar = False, title = title, vmin = np.nanmin(plot_variable_2d), vmax = np.nanmax(plot_variable_2d), cmap = plt.cm.get_cmap('plasma'), path = filename)
            #PLOT 1: rectangular plot of field open vs. closed
            title = r'{} Field line topology, time={}'.format(self.run.upper(), self.fileIndex)
            filename = '{}rect_frame_{}_fieldopen{}.png'.format(save_dir, str(self.fileIndex).zfill(5), self.file_suffix)
            self.make_lonlat_plot(phi, theta, plot_variable_2d, xlabel='longitude [deg.]', ylabel='latitude [deg.]', title = title, cbar = False,
                                  #cbar_label='Open [1] vs. closed regions [0]', 
                                  extent = [-180, 180, -90, 90],
                                  vmin = np.nanmin(plot_variable_2d), vmax = np.nanmax(plot_variable_2d), cmap = plt.cm.get_cmap('plasma'), path = filename)
            #PLOT 2: look at field lines near the noon longitudinal line
            # calculate last open and last closed lines at noon
            lon_1d = phi[0,:]
            ind = np.where(np.abs(lon_1d) == np.min(np.abs(lon_1d)))[0][0]         # find the noon meridian            
            last_open_index_list, last_closed_index_list = open_closed_boundary( plot_variable_2d[:, ind] )
            coord_list_open = [self.coord_list[0][last_open_index_list, ind], self.coord_list[1][last_open_index_list, ind], self.coord_list[2][last_open_index_list, ind] ]
            coord_list_closed = [self.coord_list[0][last_closed_index_list, ind], self.coord_list[1][last_closed_index_list, ind], self.coord_list[2][last_closed_index_list, ind] ]
            lines_open = static_field_tracer_3d_alt( self.f, coord_list_open, self.max_iterations, self.dx, direction='+-', fg_b = self.fg_b )   #list of [nlines, 3] position arrays
            lines_closed = static_field_tracer_3d_alt(self.f, coord_list_closed, self.max_iterations, self.dx, direction='+-', fg_b = self.fg_b )
            # plot projection of open & closed lines in x-z plane
            for i in range(len(last_open_index_list)):
                plt.plot(np.array(lines_open)[:, i, 0] / R_EARTH, np.array(lines_open)[:, i, 2] / R_EARTH, color = 'green')
            for j in range(len(last_closed_index_list)):
                plt.plot(np.array(lines_closed)[:, j, 0] / R_EARTH, np.array(lines_closed)[:, j, 2] / R_EARTH, color = 'red')
            # plot B field vectors
            yeq0_ind = int(fg_b.shape[1]/2)
            B = (fg_b[:,:,:,0]**2 + fg_b[:,:,:,1]**2 + fg_b[:,:,:,2]**2)**0.5
            x, y, z = fg_grid(self.f, fg_b = self.fg_b)
            x2_temp, z2_temp = np.meshgrid(x, z, indexing='ij', sparse=True)
            x2d = x2_temp + (z2_temp * 0)
            z2d = z2_temp + (x2_temp * 0)
            dq = 40        # index spacing between quivers
            plt.quiver(x2d[0::dq,0::dq]/R_EARTH, z2d[0::dq,0::dq]/R_EARTH, 
                       2*self.fg_b[0::dq,yeq0_ind,0::dq,0]/B[0::dq,yeq0_ind,0::dq], 2*self.fg_b[0::dq,yeq0_ind,0::dq,2]/B[0::dq,yeq0_ind,0::dq] )
            # plot shue model
            step = np.pi / 200   
            theta_shue = np.arange( -(np.pi-step), np.pi-step, step)
            r_shue, r_0_shue, alpha_shue = f_shue(theta_shue, run = self.run)
            plt.plot( r_shue * np.cos(theta_shue), r_shue * np.sin(theta_shue), color = 'black' )
            plt.plot(  np.cos(theta_shue), np.sin(theta_shue), color = 'blue', linestyle = 'dotted' )        # Earth has radius = 1 on this plot
            plt.title(title)
            plt.xlabel(r'x [$R_E$]')
            plt.ylabel(r'z [$R_E$]')
            plt.xlim([f.read_parameter('xmin') / R_EARTH, f.read_parameter('xmax') / R_EARTH])
            plt.ylim([f.read_parameter('zmin') / R_EARTH, f.read_parameter('zmax') / R_EARTH])
            save_dir = '{}{}/field_open/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))
            path = '{}field_open_closed_lines_{}.png'.format(save_dir, str(self.fileIndex).zfill(5))
            mkdir_path(path)
            plt.savefig(path)
            plt.close()
    


class FAC_Obj(ParamObj):
    def __init__(self, run, f, coord_list, fg_b = None, **kwargs):
        if (fg_b is None):
            fg_b = f.read_variable('fg_b')     
        self.fg_b = fg_b
        #approach 1: calculate the current all in one go
        #self.j = (1 / mu_0) * numcurl3d(fg_b, CELLSIZE_XYZ)
        self.j = fg_b*0
        #approach 2: break the calculation into bits to take up less memory (probably a bit slower)
        nx2 = int(fg_b.shape[0]/2)
        ny2 = int(fg_b.shape[1]/2)
        nz2 = int(fg_b.shape[2]/2)
        self.j[0:nx2,0:ny2,0:nz2,:] = (1 / mu_0) * numcurl3d(fg_b[0:nx2,0:ny2,0:nz2,:], CELLSIZE_XYZ)
        self.j[nx2:,0:ny2,0:nz2,:] = (1 / mu_0) * numcurl3d(fg_b[nx2:,0:ny2,0:nz2,:], CELLSIZE_XYZ)
        self.j[0:nx2,ny2:,0:nz2,:] = (1 / mu_0) * numcurl3d(fg_b[0:nx2,ny2:,0:nz2,:], CELLSIZE_XYZ)
        self.j[0:nx2,0:ny2,nz2:,:] = (1 / mu_0) * numcurl3d(fg_b[0:nx2,0:ny2,nz2:,:], CELLSIZE_XYZ)
        self.j[0:nx2,ny2:,nz2:,:] = (1 / mu_0) * numcurl3d(fg_b[0:nx2,ny2:,nz2:,:], CELLSIZE_XYZ)
        self.j[nx2:,0:ny2,nz2:,:] = (1 / mu_0) * numcurl3d(fg_b[nx2:,0:ny2,nz2:,:], CELLSIZE_XYZ)
        self.j[nx2:,ny2:,0:nz2,:] = (1 / mu_0) * numcurl3d(fg_b[nx2:,ny2:,0:nz2,:], CELLSIZE_XYZ)
        self.j[nx2:,ny2:,nz2:,:] = (1 / mu_0) * numcurl3d(fg_b[nx2:,ny2:,nz2:,:], CELLSIZE_XYZ)
        # ^^ this computation of j is SLOW and takes a TON of memory: either pre-save the results somewhere or use Markku's sidecar files
        # Use Markku's sidecar (pre-computed) files: e.g.
        # f_sidecar = pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/visualizations_2/ballooning/jlsidecar_bulk1.egl.0001758.vlsv") 
        # f_sidecar.vg_J 
        # (TODO: implement this. Need to convert getween vg and fg grids, too)
        ParamObj.__init__(self, run, f, coord_list, **kwargs) 
    def make_plot_data(self, plot_data=False, style = 'j_rad', b_scaled = True): 
        interpolators_J = interpolator_list_3d(self.f, self.j, bounds_error = False, fill_value = np.nan)
        interpolators_B = interpolator_list_3d(self.f, self.fg_b, bounds_error = False, fill_value = np.nan)
        N = self.coord_list_traced[0].size
        x_out = deepcopy(self.coord_list_traced[0]); y_out = deepcopy(self.coord_list_traced[1]); z_out = deepcopy(self.coord_list_traced[2])
        x_out = x_out.reshape(N); y_out = y_out.reshape(N); z_out = z_out.reshape(N)          # *_out are 1d arrays
        point = np.array([x_out, y_out, z_out]).T.reshape([N, 3])
        x_out_0 = deepcopy(self.coord_list[0]); y_out_0 = deepcopy(self.coord_list[1]); z_out_0 = deepcopy(self.coord_list[2])
        x_out_0 = x_out_0.reshape(N); y_out_0 = y_out_0.reshape(N); z_out_0 = z_out_0.reshape(N)          # *_out are 1d arrays
        point_0 = np.array([x_out_0, y_out_0, z_out_0]).T.reshape([N, 3])
        mask = ~(np.isnan(point[:,0] * point[:,1] * point[:,2] * point_0[:,0] * point_0[:,1] * point_0[:,2]))    
        B_x_0 = interpolators_B[0](point_0[mask,:])
        B_y_0 = interpolators_B[1](point_0[mask,:])
        B_z_0 = interpolators_B[2](point_0[mask,:])
        B_mag_0 = (B_x_0**2 + B_y_0**2 + B_z_0**2 )**0.5
        B_x = interpolators_B[0](point[mask,:])
        B_y = interpolators_B[1](point[mask,:])
        B_z = interpolators_B[2](point[mask,:])
        B_mag = (B_x**2 + B_y**2 + B_z**2 )**0.5
        B_unit_x = B_x / B_mag
        B_unit_y = B_y / B_mag
        B_unit_z = B_z / B_mag
        J_x = interpolators_J[0](point[mask,:])
        J_y = interpolators_J[1](point[mask,:])
        J_z = interpolators_J[2](point[mask,:])
        plot_variable_2d = np.zeros(N)
        if b_scaled == True:
            scaled_string = 'scaled'
            cbar_string = '(B_0 / B)'
            scaling_factor = (B_mag_0 / B_mag)
        else:
            scaled_string = ''
            cbar_string = ''
            scaling_factor = 1             # set B_mag_0 / B_mag = 0 so no scaling is performed
        if style == 'jrad':
            r_out = (x_out**2 + y_out**2 + z_out**2)**0.5
            r_unit_x = x_out / r_out ; r_unit_y = y_out / r_out ; r_unit_z = z_out / r_out
            plot_variable_2d[mask] = ((J_x * r_unit_x[mask]) + (J_y * r_unit_y[mask]) + (J_z * r_unit_z[mask])) * 1e6 * scaling_factor      #radial J, scaled by ratio of magnetic fields, i.e. div J = 0
            cbar_label = r'${}J \cdot \hat r$  $[A/km^2]$'.format(cbar_string)
            title = r'{} {} $J_r$, time={} s'.format(self.run.upper(), scaled_string, self.fileIndex)
        elif style == 'jpar':
            plot_variable_2d[mask] = ((J_x * B_unit_x) + (J_y * B_unit_y) + (J_z * B_unit_z)) * 1e6 *  scaling_factor    # parallel J, scaled by ratio of magnetic fields, assuming div J = 0
            cbar_label = r'${}J_\parallel$  $[A/km^2]$'.format(cbar_string)
            title = r'{} {} FAC, time={} s'.format(self.run.upper(), scaled_string, self.fileIndex)
        plot_variable_2d = plot_variable_2d.reshape(self.shape)
        self.data_list.append(plot_variable_2d)
        self.data_label_list.append(cbar_label)
        if plot_data:
        #PLOT 1: rectangular coordinate system
            save_dir = '{}{}/fac/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))         #CUSTOM
            filename = '{}rect_frame_{}_{}{}{}.png'.format(save_dir, str(self.fileIndex).zfill(5), style, scaled_string, self.file_suffix)
            r, theta, phi = self.coord_list_sphere
            self.make_lonlat_plot(phi, theta, plot_variable_2d, xlabel='longitude [deg.]', ylabel='latitude [deg.]', title = title, cbar_label=cbar_label, extent = [-180, 180, -90, 90],
                                  #vmin = -1.5e-8, vmax= 1.5e-8, cmap = plt.cm.get_cmap('plasma'), path = filename)
                                  vmin = np.nanmin(plot_variable_2d), vmax = np.nanmax(plot_variable_2d), cmap = 'bwr', path = filename)
            #PLOT 2: polar plot
            delta_lat = 40                   # range of latitudes (in degrees) to plot
            nlat = self.shape[0]
            delta_lat_ind = int(nlat * (delta_lat / 180))
            save_dir = '{}{}/fac/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))
            filename = '{}azim_frame_{}_{}{}{}.png'.format(save_dir, str(self.fileIndex).zfill(5), style, scaled_string, self.file_suffix)
            #mm = np.nanmax(np.abs([plot_variable_2d[(nlat-delta_lat_ind-1):,:]]))
            self.make_polar_plot(phi[(nlat-delta_lat_ind-1):,:], theta[(nlat-delta_lat_ind-1):,:], plot_variable_2d[(nlat-delta_lat_ind-1):,:],
                                 title = title, cbar_label=cbar_label, vmin = -1, vmax = 1,   # vmin = -mm, vmax = mm,
                                 #vmin = np.nanmin(plot_variable_2d[(nlat-delta_lat_ind-1):,:]), vmax = np.nanmax(plot_variable_2d[(nlat-delta_lat_ind-1):,:]),
                                 cmap = 'bwr', path = filename)
    


class MagnetopauseShockObj(ParamObj):
    '''calculate position of the magnetopause and the bow shock.
    Use that information to define the magnetosheath
    '''
    def __init__(self, run, f, coord_list, fg_b = None, **kwargs):
        if (fg_b is None):
            fg_b = f.read_variable('fg_b')     
        self.fg_b = fg_b 
        ParamObj.__init__(self, run, f, coord_list, **kwargs)
    def make_plot_data(self, plot_data=False, ignore = 3, count = 3):
        import magnetopause3dk          # don't want to do this at the top of the file because it kind of breaks vorna
        streams = magnetopause3dk.make_streamlines(self.run, self.fileIndex)
        pause_by_slice, pause_by_x = magnetopause3dk.get_magnetopause(streams, ignore, count)
        #mpos_array = np.array(streams.streamlines)    # shape (nstreams, len of stream, 3)
        mpos_array = np.array(pause_by_slice)          # shape (nslices, nx, 3) ?
        mpos_x = mpos_array[:,:,0].flatten() 
        mpos_y = mpos_array[:,:,1].flatten() 
        mpos_z = mpos_array[:,:,2].flatten() 
        mpos_r = (mpos_x**2 + mpos_y**2 + mpos_z**2)**0.5
        mpos_rcyl = (mpos_y**2 + mpos_z**2)**0.5
        theta_polar = np.arccos(mpos_x / mpos_r)
        r_0, alpha = f_shue_parameters(self.run)       # initial parameters
        params, params_cov = scipy.optimize.curve_fit(f_shue_parametrized, theta_polar, mpos_r, p0=[r_0, alpha])
        theta_polar_data = np.arccos(self.coord_list[0] / self.coord_list_sphere[0])
        plot_variable_2d = f_shue_parametrized(theta_polar_data, params[0], params[1])   # r
        plot_variable_2d = plot_variable_2d.reshape(self.shape)
        self.data_list.append(plot_variable_2d)
        self.data_label_list.append('Magnetopause_radius_RE')
        save_dir = '{}{}/magnetopause_pos/{}/'.format(self.root_dir, self.run.upper(), str(self.fileIndex).zfill(7))         #CUSTOM
        if plot_data:
            #PLOT 0: 3d mesh plot of magnetopause
            #magnetopause3dk.main(run = self.run, fileIndex = self.fileIndex, save_dir = save_dir)
            #plt.close()
            #PLOT 1: x vs r_cyl of magnetopause points
            title = r'{} Magnetopause position, time={}'.format(self.run.upper(), self.fileIndex)
            plt.scatter(mpos_x, mpos_rcyl)
            plt.xlabel(r'x [$r_E$]')
            plt.ylabel(r'$\sqrt{y^2 + z^2}$ [$r_E$]')
            plt.title(title)
            filename = '{}magnetopause_pos_data_{}_{}.png'.format(save_dir, str(self.fileIndex).zfill(5), self.file_suffix)
            mkdir_path(filename)
            plt.savefig(filename)
            plt.close()
            #PLOT 2: The fit results
            title = r'{} Shue fit $r(\theta)$, ($r_0$={}), time={}'.format(self.run.upper(), str(params[0])[0:3],  self.fileIndex)
            step = np.pi / 100
            #theta_polar_plot = np.arange(step, np.pi-step, step)
            theta_polar_plot = np.arange(-np.pi + step, np.pi-step, step)
            r_fit_mp = f_shue_parametrized(theta_polar_plot, params[0], params[1])
            plt.scatter(theta_polar, mpos_r, color = 'green')
            plt.scatter(theta_polar_plot, f_shue_parametrized(theta_polar_plot, r_0, alpha), color = 'blue')
            plt.scatter(theta_polar_plot, r_fit_mp, color='red')
            plt.xlabel(r'$\theta$')
            plt.ylabel(r'Magnetopause position [$r_E$]')
            filename = '{}rect_frame_shue_fit_{}_{}.png'.format(save_dir, str(self.fileIndex).zfill(5), self.file_suffix)
            mkdir_path(filename)
            plt.title(title)
            plt.savefig(filename)
            plt.close()
            # PLOT 3: Fit to Magnetopause using beta_prime:
            f_shue_fit_mp_beta_prime, theta_mp, r_mp = fit_magnetopause(self.f, run = self.run, root_dir = self.root_dir, fileIndex = self.fileIndex, threshold = 1, delta = 0.05, plot = True)
            r_fit_mp_beta_prime = f_shue_fit_mp_beta_prime(theta_polar_plot)
            # PLOT 4:
            title = r'{} Magnetopause and Bow shock, time={}'.format(self.run.upper(), self.fileIndex)
            #plt.scatter(mpos_r * np.cos(theta_polar), mpos_r*np.sin(theta_polar), color = 'red')
            #plt.scatter(r_mp * np.cos(theta_mp), r_mp*np.sin(theta_mp), color = 'red')
            f_shue_fit_bs_n, theta_bs_n, r_bs_n = fit_bow_shock(self.f, run = self.run, root_dir = self.root_dir, fileIndex = self.fileIndex, threshold = 2, delta = 0.05, method = 'n', plot = True)
            f_shue_fit_bs_b, theta_bs_b, r_bs_b = fit_bow_shock(self.f, run = self.run, root_dir = self.root_dir, fileIndex = self.fileIndex, threshold = 0.96, delta = 0.01, method = 'b', plot = True)
            r_fit_bs_n = f_shue_fit_bs_n(theta_polar_plot)
            r_fit_bs_b = f_shue_fit_bs_b(theta_polar_plot)
            #plt.scatter(r_bs * np.cos(theta_bs), r_bs * np.sin(theta_bs), color = 'blue', label = 'Bow shock')
            plt.plot( r_fit_mp * np.cos(theta_polar_plot), r_fit_mp * np.sin(theta_polar_plot), color='red', label = 'Magnetopause (streamline)' )
            plt.plot( r_fit_mp_beta_prime * np.cos(theta_polar_plot), r_fit_mp_beta_prime * np.sin(theta_polar_plot), color = 'green', label = r'Magnetopause ($\beta^\prime$)' )
            plt.plot( r_fit_bs_n * np.cos(theta_polar_plot), r_fit_bs_n * np.sin(theta_polar_plot), color = 'blue', label = r'Bow shock ($M_S$)' )
            plt.plot( r_fit_bs_b * np.cos(theta_polar_plot), r_fit_bs_b * np.sin(theta_polar_plot), color = 'black', label = r'Bow shock (B deflection)' )
            # plot B field vectors
            yeq0_ind = int(self.fg_b.shape[1]/2)
            B = (self.fg_b[:,:,:,0]**2 + self.fg_b[:,:,:,1]**2 + self.fg_b[:,:,:,2]**2)**0.5
            x, y, z = fg_grid(self.f, fg_b = self.fg_b)
            x2_temp, z2_temp = np.meshgrid(x, z, indexing='ij', sparse=True)
            x2d = x2_temp + (z2_temp * 0)
            z2d = z2_temp + (x2_temp * 0)
            if self.run == 'EGL':
                dq = 40        # index spacing between quivers
            else:
                dq = 10        # index spacing between quivers
            plt.quiver(x2d[0::dq,0::dq]/R_EARTH, z2d[0::dq,0::dq]/R_EARTH, 
                       2*self.fg_b[0::dq,yeq0_ind,0::dq,0]/B[0::dq,yeq0_ind,0::dq], 2*self.fg_b[0::dq,yeq0_ind,0::dq,2]/B[0::dq,yeq0_ind,0::dq] )
            plt.xlabel(r'x [$r_E$]')
            #plt.ylabel(r'$\sqrt{y^2 + z^2}$ [$r_E$]')
            plt.ylabel(r'z [$r_E$]')
            plt.xlim( [f.read_parameter('xmin')/R_EARTH, f.read_parameter('xmax')/R_EARTH] )
            plt.ylim( [f.read_parameter('zmin')/R_EARTH, f.read_parameter('zmax')/R_EARTH] )
            plt.title(title)
            plt.legend()
            filename = '{}magnetopause_bow_shock_fits_xz_{}.png'.format(save_dir, str(self.fileIndex).zfill(5))
            plt.savefig(filename)
            plt.close()
            # PLOT 5:
            title = r'{} Magnetopause and Bow shock, time={}'.format(self.run.upper(), self.fileIndex)
            plt.xlabel(r'x [$r_E$]')
            plt.ylabel(r'y [$r_E$]')
            plt.title(title)
            plt.plot( r_fit_mp * np.cos(theta_polar_plot), r_fit_mp * np.sin(theta_polar_plot), color='red', label = 'Magnetopause (streamline)' )
            plt.plot( r_fit_mp_beta_prime * np.cos(theta_polar_plot), r_fit_mp_beta_prime * np.sin(theta_polar_plot), color = 'green', label = r'Magnetopause ($\beta^\prime$)' )
            plt.plot( r_fit_bs_n * np.cos(theta_polar_plot), r_fit_bs_n * np.sin(theta_polar_plot), color = 'blue', label = r'Bow shock (density)' )
            plt.plot( r_fit_bs_b * np.cos(theta_polar_plot), r_fit_bs_b * np.sin(theta_polar_plot), color = 'black', label = r'Bow shock (B deflection)' )
            # plot B field vectors
            zeq0_ind = int(self.fg_b.shape[2]/2)
            x2_temp, y2_temp = np.meshgrid(x, y, indexing='ij', sparse=True)
            x2d = x2_temp + (y2_temp * 0)
            y2d = y2_temp + (x2_temp * 0)
            plt.quiver(x2d[0::dq,0::dq]/R_EARTH, y2d[0::dq,0::dq]/R_EARTH, 
                       2*self.fg_b[0::dq,0::dq,zeq0_ind,0]/B[0::dq,0::dq,zeq0_ind], 2*self.fg_b[0::dq,0::dq,zeq0_ind,1]/B[0::dq,0::dq,zeq0_ind] )
            plt.legend()
            plt.xlim( [f.read_parameter('xmin')/R_EARTH, f.read_parameter('xmax')/R_EARTH] )
            plt.ylim( [f.read_parameter('ymin')/R_EARTH, f.read_parameter('ymax')/R_EARTH] )
            filename = '{}magnetopause_bow_shock_fits_xy_{}.png'.format(save_dir, str(self.fileIndex).zfill(5))
            plt.savefig(filename)
            plt.close()
            # PLOT 6:  Final inferred values with Shue model
            r_magnetopause_standoff = f_shue_fit_mp_beta_prime(0)
            title = r'{} Magnetopause ($R_M$={}) and Bow shock, time={}'.format(self.run.upper(), r_magnetopause_standoff, self.fileIndex)
            #plt.scatter(mpos_r * np.cos(theta_polar), mpos_r*np.sin(theta_polar), color = 'red')
            #plt.scatter(r_mp * np.cos(theta_mp), r_mp*np.sin(theta_mp), color = 'red')
            f_shue_fit_bs_n, theta_bs_n, r_bs_n = fit_bow_shock(self.f, run = self.run, root_dir = self.root_dir, fileIndex = self.fileIndex, threshold = 2, delta = 0.05, method = 'n', plot = True)
            f_shue_fit_bs_b, theta_bs_b, r_bs_b = fit_bow_shock(self.f, run = self.run, root_dir = self.root_dir, fileIndex = self.fileIndex, threshold = 0.96, delta = 0.01, method = 'b', plot = True)
            r_fit_bs_n = f_shue_fit_bs_n(theta_polar_plot)
            r_fit_bs_b = f_shue_fit_bs_b(theta_polar_plot)
            #plt.scatter(r_bs * np.cos(theta_bs), r_bs * np.sin(theta_bs), color = 'blue', label = 'Bow shock')
            #plt.plot( r_fit_mp * np.cos(theta_polar_plot), r_fit_mp * np.sin(theta_polar_plot), color='red', label = 'Magnetopause (streamline)' )
            plt.plot( r_fit_mp_beta_prime * np.cos(theta_polar_plot), r_fit_mp_beta_prime * np.sin(theta_polar_plot), color = 'green', label = r'Magnetopause ($\beta^\prime$)' )
            plt.plot( r_fit_bs_n * np.cos(theta_polar_plot), r_fit_bs_n * np.sin(theta_polar_plot), color = 'blue', label = r'Bow shock (density)' )
            #plt.plot( r_fit_bs_b * np.cos(theta_polar_plot), r_fit_bs_b * np.sin(theta_polar_plot), color = 'black', label = r'Bow shock (B deflection)' )
            # Shue model
            step = np.pi / 200
            theta_shue = np.arange( -(np.pi-step), np.pi-step, step)
            r_shue, r_0_shue, alpha_shue = f_shue(theta_shue, run = self.run)
            plt.plot( r_shue * np.cos(theta_shue), r_shue * np.sin(theta_shue), color = 'black', label = 'Shue (1997)' )
            # plot B field vectors
            yeq0_ind = int(self.fg_b.shape[1]/2)
            B = (self.fg_b[:,:,:,0]**2 + self.fg_b[:,:,:,1]**2 + self.fg_b[:,:,:,2]**2)**0.5
            x, y, z = fg_grid(self.f, fg_b = self.fg_b)
            x2_temp, z2_temp = np.meshgrid(x, z, indexing='ij', sparse=True)
            x2d = x2_temp + (z2_temp * 0)
            z2d = z2_temp + (x2_temp * 0)
            if self.run == 'EGL' or self.run == 'EGP':
                dq = 40        # index spacing between quivers
            else:
                dq = 10        # index spacing between quivers
            plt.quiver(x2d[0::dq,0::dq]/R_EARTH, z2d[0::dq,0::dq]/R_EARTH, 
                       2*self.fg_b[0::dq,yeq0_ind,0::dq,0]/B[0::dq,yeq0_ind,0::dq], 2*self.fg_b[0::dq,yeq0_ind,0::dq,2]/B[0::dq,yeq0_ind,0::dq] )
            plt.xlabel(r'x [$r_E$]')
            #plt.ylabel(r'$\sqrt{y^2 + z^2}$ [$r_E$]')
            plt.ylabel(r'z [$r_E$]')
            plt.xlim( [-80, 40] )
            plt.ylim( [-60, 60] )
            #plt.xlim( [f.read_parameter('xmin')/R_EARTH, f.read_parameter('xmax')/R_EARTH] )
            #plt.ylim( [f.read_parameter('zmin')/R_EARTH, f.read_parameter('zmax')/R_EARTH] )
            plt.title(title)
            plt.legend()
            filename = '{}magnetopause_bow_shock_shue_xy_{}.png'.format(save_dir, str(self.fileIndex).zfill(5))
            plt.savefig(filename)
            plt.close()
            #PLOT 7: rectangular coordinate system
            title = r'{} Magnetopause Radius [$r_E$], time={}'.format(self.run.upper(), self.fileIndex)
            cbar_label =  r'{} Magnetopause Radius [$r_E$]'
            filename = '{}rect_frame_magnetopause_pos_{}_{}.png'.format(save_dir, str(self.fileIndex).zfill(5), self.file_suffix)
            r, theta, phi = self.coord_list_sphere
            self.make_lonlat_plot(phi, theta, plot_variable_2d, xlabel='longitude [deg.]', ylabel='latitude [deg.]', title = title, cbar_label=cbar_label, extent = [-180, 180, -90, 90],
                                  vmin = 0, vmax = 100, path = filename)
            


#########################################################
#main routine, pass command line arguments as a dictionary
#########################################################

def analyze_data(fileIndex_list):

    if type(fileIndex_list) == int:
        fileIndex_list = [fileIndex_list]

    print(fileIndex_list)

    for fileIndex in fileIndex_list:
    
        # load data

        filename = get_vlsvfile_fullpath(run, fileIndex)
        f = pt.vlsvfile.VlsvReader( filename )
        r_aurora = 1. * R_EARTH
        r_trace = 5 * R_EARTH     
        dx = R_EARTH / 50
        max_iterations = int(4 * r_trace / dx)
        
        # make longitude-latitude grid over which to evaluate data
        #nlat = 180
        #nphi = 360
        #lat, phi = lat_phi_grid( nlat = nlat, nphi = nphi)
        degtorad = np.pi / 180
        lat, phi = lat_phi_grid( phi_min = float(ARGS.phimin)*degtorad, phi_max = float(ARGS.phimax)*degtorad, 
                                 lat_min = float(ARGS.latmin)*degtorad, lat_max = float(ARGS.latmax)*degtorad,
                                 nlat = int(ARGS.nlat), nphi = int(ARGS.nphi))
        theta = lat2theta(lat)
        coord_list_aurora = [(phi*0) + r_aurora, theta, phi]
        coord_list_aurora_cart = spherical_to_cartesian( *coord_list_aurora )
        coord_list_magnetosphere = [(phi*0) + r_trace, theta, phi] 
        coord_list_magnetosphere_cart = spherical_to_cartesian( *coord_list_magnetosphere )
        
        if ('1' in varnames) or ('2' in varnames) or ('3' in varnames) or ('4' in varnames) or ('5' in varnames):
        
            # analysator
            fg_b = f.read_variable('fg_b')     # this is SLOW  ... TODO: implement hongyong's julia version  
        
            # julia (reads fg_b faster)
            #from juliacall import Main as jl
            #jl.seval("using Vlasiator")
            #meta = jl.load(filename)
            #fg_b = jl.readvariable(meta, "fg_b")
            #fg_b = np.transpose(fg_b, (1,2,3,0))
        
            global CELLSIZE_XYZ
            #CELLSIZE = (f.read_parameter('xmax') - f.read_parameter('xmin')) / fg_b.shape[0]
            CELLSIZE_XYZ = [ (f.read_parameter('xmax') - f.read_parameter('xmin')) / fg_b.shape[0], 
                            (f.read_parameter('ymax') - f.read_parameter('ymin')) / fg_b.shape[1], 
                            (f.read_parameter('zmax') - f.read_parameter('zmin')) / fg_b.shape[2] ]
        
            # trace lines in a radial orientation
            trace_method = 'integrateB'            #options: 'dipole', 'integrateB'
            coord_list_traced_pos = trace_coord(f, coord_list_aurora, fg_b = fg_b, trace_method = trace_method, direction = '+', 
                                                coord_in = 'spherical', coord_out = 'spherical', r_trace = r_trace, max_iterations = max_iterations, dx = dx)
            coord_list_traced_neg = trace_coord(f, coord_list_aurora, fg_b = fg_b, trace_method = trace_method, direction = '-',
                                                coord_in = 'spherical', coord_out = 'spherical', r_trace = r_trace, max_iterations = max_iterations, dx = dx)
            mask = test_radial_field(f, coord_list_aurora_cart)
            coord_list_traced = coord_list_traced_neg
            for i in range(3):
                coord_list_traced[i][mask] = coord_list_traced_pos[i][mask]
        
                
            save_list = []
            # make plots and save the data for different parameters in objects 
            if '1' in varnames:
                # 1. open vs. closed field boundary
                fOO = FieldOpenObj(run, f, coord_list_aurora, fileIndex = fileIndex, fg_b = fg_b, max_iterations = 4000, dx = R_EARTH / 10, make_plot_data = True, plot_data=ARGS.plot)
                save_list.append(fOO)
            if '2' in varnames:
                # 2. proton differential energy flux
                #pFO = pFluxObj(run, f, coord_list_magnetosphere, make_plot_data = True )
                if run != 'EGP':     # EGP doesn't have particle flux data
                    pFO_traced = pFluxObj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, make_plot_data = True, plot_data=ARGS.plot)
                save_list.append(pFO_traced)
            if '3' in varnames:
                # 3. Field aligned currents (FACS)
                #JO = FAC_Obj(run, f, coord_list_magnetosphere, fileIndex = fileIndex, fg_b = fg_b)
                #JO.make_plot_data(style = 'jpar')
                JO_traced = FAC_Obj(run, f, coord_list_aurora, fileIndex = fileIndex, coord_list_traced = coord_list_traced, fg_b = fg_b, plot_data=ARGS.plot)
                JO_traced.make_plot_data(style = 'jpar', b_scaled = True)
                save_list.append(JO_traced)
            if '4' in varnames:
                # 4. Magnetopause position
                MSO = MagnetopauseShockObj(run, f, coord_list_aurora, fileIndex = fileIndex, fg_b=fg_b, make_plot_data = True, plot_data=ARGS.plot)
                save_list.append(MSO)
                # merge objects and write data to a .csv file, to be shared with FMI collaborators
                # note: all data evaluated at the same coordinates
            if '5' in varnames:
                # 5. dB/dt. Note: this needs to be evaluated in the ionosphere (once fully implemented by Urs)
                dBO = dBdtObj(run, f, coord_list_magnetosphere, fileIndex = fileIndex,  make_plot_data = True, file_suffix = '_5re', plot_data=ARGS.plot)
                save_list.append(dBO)
            if ARGS.save == True:
                # Save data into a .csv file
                save_obj = ParamObj(run, f, coord_list_aurora, coord_list_traced = coord_list_traced)
                save_obj.merge(save_list)
                save_obj.write_data( filename = ARGS.savefile, append=ARGS.append)
        
        # 6. other plots
        if '6' in varnames:
            save_dir = '{}{}/summary/{}/'.format(ROOT_DIR, run.upper(), str(fileIndex).zfill(7))         #CUSTOM
            mkdir_path(save_dir)
    #            pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_rho',normal='y', boxre=[-30,30,-30,30], vmin = 1e4, vmax = 1e8, run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_rho_yeq0_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_b_vol', cutpointre=0)
    #            pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_v',operator='x',normal='y', boxre=[-40,0,-15,15], vmin = -2e6, vmax = 2e6, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_vx_yeq0_tail_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_b_vol',streamlinecolor='black',streamlinedensity=2,streamlinethick=0.8, cutpointre=0)
    #            pt.plot.plot_colormap3dslice(filename=filename,var='vg_b_vol',operator='z',normal='y', boxre=[-40,0,-15,15], vmin = -1.5e-8, vmax = 1.5e-8, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_bz_yeq0_tail_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_b_vol',streamlinecolor='black',streamlinedensity=2,streamlinethick=0.8, cutpointre=0)
    #            pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_rho',normal='y', boxre=[-15,15,-15,15],  vmin = 1e5, vmax = 1e8,run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_rho_gradpe_lines_yeq0_tail_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_e_gradpe',streamlinecolor='black',streamlinedensity=2,streamlinethick=0.8, cutpointre=0)
            if run != 'EGI':
                pt.plot.plot_colormap3dslice(filename=filename,var='vg_e_parallel',normal='y', boxre=[-20,0,-15,15], vmin = -2e-3, vmax = 2e-3, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_e_parallel_yeq0_tail_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_b_vol',streamlinecolor='black',streamlinedensity=2,streamlinethick=0.8, cutpointre=0)
    #                pt.plot.plot_colormap3dslice(filename=filename,var='vg_poynting',normal='y', boxre=[-20,0,-15,15], vmin = 1e-5, vmax = 1e-2, run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_poynting_yeq0_tail_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_b_vol',streamlinecolor='black',streamlinedensity=2,streamlinethick=0.8, cutpointre=0)
    #                pt.plot.plot_colormap3dslice(filename=filename,var='vg_poynting',operator='x',normal='y', boxre=[-20,0,-15,15], vmin = -5e-3, vmax = 5e-3, run=run,colormap='bwr',step=fileIndex,outputdir=save_dir,outputfile='vg_poynting_x_yeq0_tail_{}_{}.png'.format(run,fileIndex),Earth=1, streamlines = 'vg_b_vol',streamlinecolor='black',streamlinedensity=2,streamlinethick=0.8, cutpointre=0)
            #for channel in range(16):
    #            for channel in range(9):
    #               proton_energy = str(int(f.read_parameter('proton_PrecipitationCentreEnergy{}'.format(channel)))).zfill(5)
    #               pt.plot.plot_colormap3dslice(filename=filename,var='proton/vg_precipitationdifferentialflux', operator = '{}'.format(channel), normal='y', boxre=[-30,30,-30,30], vmin = 1e0, vmax = 1e6, run=run,colormap='plasma',step=fileIndex,outputdir=save_dir,outputfile='vg_precipitationdifferentialflux_yeq0_{}ev_{}_{}.png'.format(proton_energy,run,fileIndex),Earth=1, streamlines = 'vg_b_vol',streamlinecolor='palegreen',streamlinethick=0.8, cutpointre=0)
        



 
if __name__ == "__main__":

    time_start = time() # measure time it takes to run the program

    # example call from bash (see carrington.sh):
    # python carrington.py -run EGI -var 1 4 6

    # Input parameters
    parser = argparse.ArgumentParser()

    parser.add_argument('-nproc', default=1, help="number of processors to use " )
    parser.add_argument('-startstop', nargs='*', help="2-element list, start and stop index (divided by deltanframes)" )
    parser.add_argument('-deltanframes',  help="only analyze one in every delta_nframes" )
    parser.add_argument('-run', default='EGL', help="the Vlasiator run, in all caps " )
    parser.add_argument('-var', default = ['1','2','3','4','5','6'], nargs='*', help="a list of plot identifiers (numbers), set to which plots you want to make" )
    parser.add_argument('-nphi', default=360, help="the number of longitudes in the lon-lat grid" )
    parser.add_argument('-phimin', default = -180, help="minimum longitude, (range -180 to 180)" ) # note phi = (pi/180) * (phi_input-180)
    parser.add_argument('-phimax', default = 180, help="maximum longitude, (range -180 to 180)" )
    parser.add_argument('-nlat', default=180, help="the number of latitudes in the lon-lat grid" )
    parser.add_argument('-latmin', default = -90, help="minimum longitude, (range -90 to 90)" ) # note lat = -lat_input * (pi/180) 
    parser.add_argument('-latmax', default = 90, help="maximum longitude, (range -90 to 90)" )
    parser.add_argument('-plot', action='store_true', help="set this flag to save the data into a .csv format")  #default: plot=False when -save flag not set
    parser.add_argument('-save', action='store_true', help="set this flag to save the data into a .csv format")  #default: save=False when -save flag not set
    parser.add_argument('-savefile', default='/wrk-vakka/users/horakons/carrington/data/test_data.csv', help="filename for the .csv data")
    parser.add_argument('-append', action='store_true', help="set this flag to append to .csv file")  #default: append=False when -save flag not set

    global ARGS
    ARGS = parser.parse_args()
    varnames = ARGS.var 

    #julia
    #if ('1' in varnames) or ('2' in varnames) or ('3' in varnames) or ('4' in varnames) or ('5' in varnames):
    #    from juliacall import Main as jl

    run = ARGS.run

    # for some runs, make one plot for only 1 every delta_nframes times
    if ARGS.deltanframes is None:
        if run == 'EGL' or run == 'EGI':
            delta_nframes = 20
        elif run == 'EGP':
            delta_nframes = 1
        else:
            delta_nframes = 1
    else:
        delta_nframes = int(ARGS.deltanframes)


#    # Frame extent for this job given as command-line arguments
    if ARGS.startstop is None:
        if run == 'EGI':
            fileIndex_list = [1500]     # 662-1506  EGI (note dB/dt is broken for t=0)
        elif run == 'EGL':
            fileIndex_list = [1760]     # EGL
        elif run == 'EGM':
            fileIndex_list = [1247]    # EGM
        elif run == 'EGN':
            fileIndex_list = [488]     # EGN
        elif run == 'EGO':
            fileIndex_list = [151]     # EGO
        elif run == 'EGP':
            fileIndex_list = [319]     # bulk1, EGP (available fileIndex: 269-319)
            #fileIndex_list = [299]     # bulk1, EGP (available fileIndex: 269-319)
            #fileIndex_list = [53]     # bulk5, EGP (available fileIndex: 1-53)
    else:
        #fileIndex_list = range(int(ARGS.startstop[0]), int(ARGS.startstop[1])+1, delta_nframes)
        fileIndex_list = range(int(ARGS.startstop[0])*delta_nframes, (int(ARGS.startstop[1]))*delta_nframes, delta_nframes)


    print(fileIndex_list)

    ## Parallel processing
    from multiprocessing import Pool
    pool = Pool(int(ARGS.nproc))
    return_array = pool.map(analyze_data, fileIndex_list)


    # measure time it takes to run the program
    time_stop = time() 
    print('total time is: {} seconds'.format(time_stop - time_start))




