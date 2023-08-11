##!/appl/easybuild/opt/Python/3.8.2-GCCcore-9.3.0/bin/python

import pytools as pt
import numpy as np
from myutils import spherical_to_cartesian, cartesian_to_spherical, get_vlsvfile_fullpath, timer, save, mkdir_path, restore
from carrington import get_all_cell_coordinates
import matplotlib.pyplot as plt
from numba import jit
from memory_profiler import profile      # @profile decorator

global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations
#global CELLSIZE
global ROOT_DIR
ROOT_DIR = '/wrk-vakka/users/horakons/carrington/plots/'
#ROOT_DIR = '/wrk-vakka/users/horakons/carrington/test_plots/'
global mu_0
mu_0 = 4e-7 * np.pi



# Input parameters
import argparse
parser = argparse.ArgumentParser()
    
parser.add_argument('-run', default='EGL', help="run name" )
parser.add_argument('-task', default=0, help="task no." )
parser.add_argument('-nproc', default=1, help="number of processors to use " )
global ARGS 
ARGS = parser.parse_args()
#nproc = ARGS.nproc



def fac_map(f, vg_x, vg_y, vg_z, f_J_sidecar = None, r_io = 5 * 6.371e6, dx_1re = 2e6, mag_mom_vector = np.array([0., 0., -8e22])):
    '''
     Find the vector current density [A/m^2] at specified x,y,z position [SI]
     Map the FACs along magnetic field lines (J \propto B).
     whether the FACs are mapped down to the ionosphere or up to the inner magnetosphere boundary depends on whether f_J_sidecar keyword is set
     f: VlsvReader object
     f_J_sidecar: vlsvReader object that contains pre-computed current 'vg_J'
     e.g., for EGL, files at: /wrk-vakka/group/spacephysics/vlasiator/3D/EGL/visualizations_2/ballooning/*.vlsv
     vg_x,vg_y,vg_z position [m], 1D numpy arrays.
     dx_1re is the grid resolution of the fac region at r= R_EARTH. This is for the input cells which can be defined arbitrarily---not necessarily same as vg cells.
     *** Works for arbitrary x,y,z coordinates, not just coordinates on the vg_ grid ***
    '''
    cellids = f.read_variable('CellID')
    dx = ((f.read_parameter('xmax') - f.read_parameter('xmin')) / f.read_parameter('xcells_ini')) # refinement lvl 0 (coarsest)
    vg_b_vol = b_dip(vg_x, vg_y, vg_z, mag_mom_vector = mag_mom_vector)

    vg_r, vg_theta, vg_phi = cartesian_to_spherical(vg_x, vg_y, vg_z)
    #ind_r = np.where(vg_r < R_EARTH)[0]
    #dx_1re = dx / 2**np.nanmin(f.get_amr_level(cellids[ind_r]))   # grid resolution dx at radius 1 RE (well inside inner boundary)
    #outer, = np.where(vg_r >= r_io)
    inner, = np.where((vg_r < r_io) & (vg_r > (R_EARTH + dx_1re/2)))  # only integrate over cells whose centers are at least 1/2 a cell length beyond radius of 1 R_E


    vg_J_eval = np.zeros([vg_x.size, 3])
    #vg_J_eval = np.zeros(vg_J.shape) # the currents that will actually be integrated over
    #vg_J_eval[outer] = vg_J[outer]

    vg_lat = (np.pi / 2) - vg_theta         # -pi/2 < lat < pi/2
    L = vg_r /  np.cos(vg_lat)**2           # L-shell (dipole) [m]

    # evaluate FACs in 'inner' region (outside of simulation domain)
    if f_J_sidecar is None:
        # map: initial point -> downmap to ionosphere via dipole formula (ionospheric runs, e.g. FHA)
        coords_ionosphere = f.get_ionosphere_node_coords()  # [n_nodes, 3]
        ig_fac = f.read_variable('ig_fac')      # [n_nodes]
        vg_b_vol_magnitude = np.sqrt(vg_b_vol[:,0]**2 + vg_b_vol[:,1]**2 + vg_b_vol[:,2]**2 )

        lat0 = np.arccos( np.sqrt(R_EARTH / L) ) # latitude at r=R_EARTH
        theta0 = (np.pi / 2) - lat0
        b0 = b_dip_magnitude(theta0, R_EARTH, mag_mom = 8e22)
        x0, y0, z0 = spherical_to_cartesian(R_EARTH, theta0[inner], vg_phi[inner])
        for i in range(x0.size):
            # find the nearest cell and evaluate the current there (brute force)
            # this approach is probably faster: https://github.com/fmihpc/vlasiator/blob/master/sysboundary/ionosphere.cpp#L381
            dist = np.sqrt((x0[i] - coords_ionosphere[:,0])**2 + (y0[i] - coords_ionosphere[:,1])**2 + (z0[i] - coords_ionosphere[:,2])**2)
            ind_min = np.argmin(dist)
            vg_J_eval[inner[i], :] = (vg_b_vol[inner[i],:] / b0[inner[i]]) * ig_fac[ind_min]  # J \propto B. Mapping UP from the FACs evaluated at the ground 
    else: # (use sidecar containing current density "vg_J" in non-ionospheric runs, e.g. EGL)
        # map: initial point -> some point in the simulation domain near the inner boundary (~5 R_E) according to dipole formula
        print('NOTE: Upmapping FACs along constant L-shell via dipole formula!')
        r_up = r_io
        lat_up = np.arccos( np.sqrt(r_up / L) ) # latitude at r=r_up
        theta_up = (np.pi / 2) - lat_up
        x_up, y_up, z_up = spherical_to_cartesian(r_up, theta_up[inner], vg_phi[inner])
        ind_fin, = np.where(np.isfinite(lat_up[inner]))
        coords_temp = list(np.array([x_up[ind_fin],y_up[ind_fin],z_up[ind_fin]]).T.reshape([ind_fin.size,3]))
        vg_b_vol_fin = f_J_sidecar.read_interpolated_variable("vg_b_vol", coords_temp)
        B_up = np.array([np.linalg.norm(vg_b_vol_fin, axis = 1)] * 3).transpose()
        B_down = np.array([b_dip_magnitude(vg_theta[inner[ind_fin]], vg_r[inner[ind_fin]], mag_mom = 8e22)] * 3).transpose()
        scale_factor = B_down / B_up                                        # J \propto B
        vg_J = f_J_sidecar.read_interpolated_variable("vg_J", coords_temp)
        J_signed_up = np.array([np.sum(vg_J * vg_b_vol_fin, axis = 1) ] * 3).transpose() / B_up    # magnitude and sign of J   (projection J dot B / |B|)
        b_dir = b_dip_direction(vg_x[inner[ind_fin]], vg_y[inner[ind_fin]], vg_z[inner[ind_fin]])
        vg_J_eval[inner[ind_fin], :] = b_dir * J_signed_up * scale_factor   # Mapping DOWN from the FACs evaluated in the simulation domain near inner boundary

    # don't allow field aligned currents if they can't possibly map to the inner boundary (some sort of numerical error led to non-zero ig_fac at equatorial latitudes)
    ind_to_zero = np.where(L < r_io) 
    vg_J_eval[ind_to_zero, :] = 0.
    return vg_J_eval


@timer
def biot_savart(coord_list, f, f_J_sidecar = None, r_io = 5 * 6.371e6):
    '''
    param coord_list:   a list of 3-element arrays of coordinates [ [x1,y1,z1], [x2,y2,z2], ... ], SI units
                        if considering just a single starting point, the code accepts a 3-element array-like object [x1,y1,z1]
    f: vlsvReader object
    f_J_sidecar: vlsvReader object that contains pre-computed current 'vg_J'
        e.g., for EGL, files at: /wrk-vakka/group/spacephysics/vlasiator/3D/EGL/visualizations_2/ballooning/*.vlsv

    runtime (FHA): overhead of about 200 sec (setup), plus 0.2 sec for each element of coord_list

    returns a tuple (B_inner, B_outer) 
    B_inner: the B-field [T] generated by FACs at 1 R_E < r< r_io
    B_outer: the B-field [T] generated by currents in simulation domain at r > r_io
    '''

    # standardize input (a list of 3-element arrays/lists)

    if type(coord_list[0]) not in [list, np.ndarray]:
        coord_list = [coord_list]

    ncoords = len(coord_list)

    # load vg_ coordinates and compute cell volumes (dV)

    vg_b_vol = f.read_variable('vg_b_vol')
    vg_x, vg_y, vg_z = get_all_cell_coordinates(f) 

    cellids = f.read_variable('CellID')
    ncells = cellids.size
    dx = ((f.read_parameter('xmax') - f.read_parameter('xmin')) / f.read_parameter('xcells_ini')) # refinement lvl 0 (coarsest)
    dV = np.zeros(ncells)
    for i, cellid in enumerate(cellids):
        dV[i] = dx**3 / 2**(3*f.get_amr_level(cellid))

    # load or calculate currents

    if f_J_sidecar is None:
        # calculate directly from B-field Jacobian (ionospheric runs, e.g. FHA)
        vg_J = np.zeros([ncells,3])
        vg_dperbxvoldy = f.read_variable('vg_dperbxvoldy')
        vg_dperbyvoldx = f.read_variable('vg_dperbyvoldx')
        vg_dperbxvoldz = f.read_variable('vg_dperbxvoldz')
        vg_dperbzvoldx = f.read_variable('vg_dperbzvoldx')
        vg_dperbyvoldz = f.read_variable('vg_dperbyvoldz')
        vg_dperbzvoldy = f.read_variable('vg_dperbzvoldy')
        vg_J[:,0] = (1 / mu_0) * (vg_dperbzvoldy - vg_dperbyvoldz)
        vg_J[:,1] = (1 / mu_0) * (vg_dperbxvoldz - vg_dperbzvoldx)
        vg_J[:,2] = (1 / mu_0) * (vg_dperbyvoldx - vg_dperbxvoldy)
        # FILL IN: if above doesn't work, calculate J from fg_b and resample to vg grid
    else:
        vg_J = f_J_sidecar.read_variable('vg_J')

    # compute B at 'coord_list' points according to Biot-Savart law (accelerate with numba?)

    #B = np.zeros([ncoords, 3])
    
    vg_r, vg_theta, vg_phi = cartesian_to_spherical(vg_x, vg_y, vg_z)
    ind_r = np.where(vg_r < R_EARTH)[0]
    dx_1re = dx / 2**np.nanmin(f.get_amr_level(cellids[ind_r]))   # grid resolution dx at radius 1 RE (well inside inner boundary)
    outer, = np.where(vg_r >= r_io)
    inner, = np.where((vg_r < r_io) & (vg_r > (R_EARTH + dx_1re/2)))  # only integrate over cells whose centers are at least 1/2 a cell length beyond radius of 1 R_E

    vg_J_eval = np.zeros(vg_J.shape) # the currents that will actually be integrated over
    vg_J_eval[outer] = vg_J[outer]

    vg_J_eval[inner] = fac_map(f, vg_x[inner], vg_y[inner], vg_z[inner], f_J_sidecar = f_J_sidecar, r_io = r_io, dx_1re = dx_1re, mag_mom_vector = np.array([0., 0., -8e22]))

    '''
    vg_lat = (np.pi / 2) - vg_theta         # -pi/2 < lat < pi/2
    L = vg_r /  np.cos(vg_lat)**2           # L-shell (dipole) [m]

    # evaluate FACs in 'inner' region (outside of simulation domain)
    if f_J_sidecar is None:
        # map: initial point -> downmap to ionosphere via dipole formula (ionospheric runs, e.g. FHA)
        coords_ionosphere = f.get_ionosphere_node_coords()  # [n_nodes, 3]
        ig_fac = f.read_variable('ig_fac')      # [n_nodes]
        vg_b_vol_magnitude = np.sqrt(vg_b_vol[:,0]**2 + vg_b_vol[:,1]**2 + vg_b_vol[:,2]**2 )

        lat0 = np.arccos( np.sqrt(R_EARTH / L) ) # latitude at r=R_EARTH
        theta0 = (np.pi / 2) - lat0
        b0 = b_dip_magnitude(theta0, R_EARTH, mag_mom = 8e22)
        x0, y0, z0 = spherical_to_cartesian(R_EARTH, theta0[inner], vg_phi[inner])
        for i in range(x0.size):
            # find the nearest cell and evaluate the current there (brute force)
            # this approach is probably faster: https://github.com/fmihpc/vlasiator/blob/master/sysboundary/ionosphere.cpp#L381
            dist = np.sqrt((x0[i] - coords_ionosphere[:,0])**2 + (y0[i] - coords_ionosphere[:,1])**2 + (z0[i] - coords_ionosphere[:,2])**2)
            ind_min = np.argmin(dist)
            vg_J_eval[inner[i], :] = (vg_b_vol[inner[i],:] / b0[inner[i]]) * ig_fac[ind_min]  # J \propto B. Mapping UP from the FACs evaluated at the ground 
    else: # (use sidecar containing current density "vg_J" in non-ionospheric runs, e.g. EGL)
        # map: initial point -> some point in the simulation domain near the inner boundary (~5 R_E) according to dipole formula
        print('NOTE: Upmapping FACs along constant L-shell via dipole formula!')
        r_up = r_io
        lat_up = np.arccos( np.sqrt(r_up / L) ) # latitude at r=r_up
        theta_up = (np.pi / 2) - lat_up
        x_up, y_up, z_up = spherical_to_cartesian(r_up, theta_up[inner], vg_phi[inner])
        ind_fin, = np.where(np.isfinite(lat_up[inner]))
        coords_temp = list(np.array([x_up[ind_fin],y_up[ind_fin],z_up[ind_fin]]).T.reshape([ind_fin.size,3]))
        vg_b_vol_fin = f_J_sidecar.read_interpolated_variable("vg_b_vol", coords_temp)
        B_up = np.array([np.linalg.norm(vg_b_vol_fin, axis = 1)] * 3).transpose()
        B_down = np.array([b_dip_magnitude(vg_theta[inner[ind_fin]], vg_r[inner[ind_fin]], mag_mom = 8e22)] * 3).transpose()
        scale_factor = B_down / B_up                                        # J \propto B
        vg_J = f_J_sidecar.read_interpolated_variable("vg_J", coords_temp)
        J_signed_up = np.array([np.sum(vg_J * vg_b_vol_fin, axis = 1) ] * 3).transpose() / B_up    # magnitude and sign of J   (projection J dot B / |B|)
        b_dir = b_dip_direction(vg_x[inner[ind_fin]], vg_y[inner[ind_fin]], vg_z[inner[ind_fin]])
        vg_J_eval[inner[ind_fin], :] = b_dir * J_signed_up * scale_factor   # Mapping DOWN from the FACs evaluated in the simulation domain near inner boundary
    '''

    # 'outer' magnetopsheric contribution given by Vlasiator currents
    # 'inner' magnetopsheric contribution given by mapped FACs

    B_inner = integrate_biot_savart(coord_list, vg_x[inner], vg_y[inner], vg_z[inner], vg_J_eval[inner], dV[inner])
    B_outer = integrate_biot_savart(coord_list, vg_x[outer], vg_y[outer], vg_z[outer], vg_J_eval[outer], dV[outer])

    return B_inner, B_outer



def b_dip_magnitude(theta, r, mag_mom = 8e22):
    # default: mag_mom = 8e22 [A / m^2] magnetic moment dipole, as in EGI, EGL, FHA runs
    B_magnitude = mag_mom * (mu_0 / (4 * np.pi * r**3)) * np.sqrt((2*np.cos(theta))**2 + (np.sin(theta))**2)
    return B_magnitude


def b_dip_direction(x, y, z, mag_mom_vector = np.array([0., 0., -8e22])):
    B = b_dip(x, y, z, mag_mom_vector = mag_mom_vector)
    return B / np.array([np.linalg.norm(B, axis = 1)] * 3).transpose()

def b_dip(x, y, z, mag_mom_vector = np.array([0., 0., -8e22])):
    N = x.size
    pos_N = np.array([x, y, z]).transpose()    # shape (N, 3)
    m_N = np.array([list(mag_mom_vector)]*N)  # shape (N, 3)
    r_N = np.array([np.linalg.norm(pos_N, axis = 1)] * 3).transpose()   # radius, shape (N, 3)
    # dipole field:  B(r) = (mu_0 / 4 pi) * (3r (m dot r) / r^5 - m / r^3)
    B = (mu_0 / (4 * np.pi)) * ( ( 3 * pos_N * np.array([np.sum(m_N * pos_N, axis = 1)]*3).transpose() / r_N**5) - m_N / r_N**3 )
    return B



@jit(nopython=True)
def integrate_biot_savart(coord_list, vg_x, vg_y, vg_z, vg_J, delta):
    # accelerated with numba
    # exact same formula can be used for volume and surface integrals
    # Biot-Savart (volume): B = (mu_0 / 4 * pi) \int { J x r' / |r'|^3 } dV  ([J] = A/m^2, delta == dV)
    #            (surface): B = (mu_0 / 4 * pi) \int { J x r' / |r'|^3 } dA  ([J] = A/m, delta = dS)
    B = np.zeros((len(coord_list), 3))
    r_p = np.zeros((vg_x.size, 3))

    for i, coord in enumerate(coord_list):
        r_p[:,0] = coord[0] - vg_x
        r_p[:,1] = coord[1] - vg_y
        r_p[:,2] = coord[2] - vg_z

        r_p_mag = np.sqrt(r_p[:,0]**2 + r_p[:,1]**2 + r_p[:,2]**2)
        J_cross_r_p = np.cross(vg_J, r_p)

        B[i,0] += np.nansum( (mu_0 / (4 * np.pi)) * delta * J_cross_r_p[:,0] / r_p_mag**3 )
        B[i,1] += np.nansum( (mu_0 / (4 * np.pi)) * delta * J_cross_r_p[:,1] / r_p_mag**3 )
        B[i,2] += np.nansum( (mu_0 / (4 * np.pi)) * delta * J_cross_r_p[:,2] / r_p_mag**3 )

    return B



#def ionosphere_node_xyz(filename = '/wrk-vakka/group/spacephysics/vlasiator/temp/ionogrid_FHA.vlsv'):
#    # load ionospheric grid (for pre-ionospheric runs, use auxiliary file)
#    # /wrk-vakka/group/spacephysics/vlasiator/temp/ionogrid_FHA.vlsv
#    f = pt.vlsvfile.VlsvReader(filename)
#    latlon = f.get_ionosphere_latlon_coords()       # 0<lat<pi, -pi<lon<pi
#    # convert to cartesian x,y,z [SI]
#    # convert to coord_list
#    x, y, z = spherical_to_cartesian(R_EARTH, latlon[:,0], latlon[:,1])
#    return list(np.array([x,y,z]).T.reshape([x.size,3]))



def get_ig_r(f):
    # ionospheric mesh is at radius (R_EARTH + 100 km), see Urs's ionosphere writeup
    n = f.get_ionosphere_node_coords()          # node = vertex of the triangular mesh
    ec = f.get_ionosphere_element_corners()     # (Element Corners), where element = trianglular face
    ig_r = np.zeros(ec.shape)
    for i in range(ig_r.shape[0]):
        ig_r[i,:] = (n[ec[i,0], :] + n[ec[i,1], :] + n[ec[i,2], :]) / 3  #barycenter, aka centroid
    return ig_r    # [n_elements, 3]



def calc_Dst(f, f_J_sidecar = None, r_io = 5 * 6.371e6):
    '''
    Follow Welling et al. (2020) approach for calculating Dst. Integrate Biot-Savart over:
    1. All currents within the Vlasov domain
    2. Birkeland currents (FACs) in the “gap region” between the MHD inner boundary and the ionosphere, 
        mapped along assumed dipole field lines
    3. Ionospheric Hall currents
    4. Ionospheric Pedersen currents
    '''
    coord_list = [np.array([R_EARTH, 0, 0]), np.array([0,R_EARTH,0]), np.array([-R_EARTH,0,0]), np.array([0,-R_EARTH,0])]  # 4 virtual magnetometers around the equator used to calculate Dst
    B_inner, B_outer = biot_savart(coord_list, f, f_J_sidecar = f_J_sidecar, r_io = r_io)    
    B = B_inner + B_outer
    Dst = np.average(B[:,2])   # for more general formula (non-equatorial), need to include cos(theta) factor
    return Dst


def ionosphere_mesh_area(f):
    # this function could be added to master
    n = f.get_ionosphere_node_coords()       # nodes: shape (21568, 3) vertices
    c = f.get_ionosphere_element_corners()   # corners of elements: indices integers 0-21567, shape (43132, 3)
    p = n[c,:]                               # shape(43132, 3, 3)   first index is the element, second index identifies corners of the triangle, third index is the x-y-z position
    r1 = p[:,1,:] - p[:,0,:] 
    r2 = p[:,2,:] - p[:,0,:] 
    areas = np.linalg.norm(np.cross(r1, r2), axis = 1) / 2.     # use cross product to calculate areas
    # checked: sum of triangle areas is near the expected area 4*pi*R^2 for a sphere
    # ( np.sum(areas) - (np.pi * 4 ) * (R_EARTH + 100000.)**2 ) / np.sum(areas) 
    return areas


def B_ionosphere(f, ig_r = None, method = 'integrate'):
    # assume currents flow around a sphere of radius 1 R_E, 
    # which locally looks like infinite plane to a ground observer looking up
    # B = (mu_0 / 2) * r_hat x J_s , where J_s vector is current per unit length
    # method = 'integrate' (default) or 'local'
    R_iono = R_EARTH + 1e5   # Ionosphere mesh at altitude 100 km, see Urs's ionosphere writeup
    if ig_r is None:
        ig_r = get_ig_r(f)
    try:
        ig_r_hat = ig_r / R_iono   # approximate (technically |ig_r| not exactly R_EARTH)
        ig_inplanecurrent = f.read_variable('ig_inplanecurrent')  # height-integrated, element centered. Units [A/m]
    except:
        return None   # no ionospheric inplanecurrent data
        #B_ionosphere = np.zeros(ig_r.shape)
    if method == "local":
        # approximate horizontal current as an infinite sheet of current directly overhead
        B_iono = (mu_0 / 2) * np.cross(ig_r_hat, ig_inplanecurrent)
    elif method == 'integrate':
        # integrate Biot-Savart law over ionospheric mesh. More accurate but slower.
        surface_r = ig_r * R_EARTH / R_iono       # same as ionosphere mesh, but at radius R_EARTH
        dS = ionosphere_mesh_area(f)
        B_iono =  integrate_biot_savart(coord_list, ig_r[:, 0], ig_r[:, 1], ig_r[:, 2], ig_inplanecurrent, dS)
    return B_iono



def B_magnetosphere(f, f_J_sidecar = None, r_io = 5 * 6.371e6, ig_r = None):
    # wrapper for biot_savart()
    if ig_r is None:
        ig_r = get_ig_r(f)
    B_inner, B_outer = biot_savart( list(ig_r), f, f_J_sidecar = f_J_sidecar, r_io = r_io )
    return B_inner, B_outer



#@profile
def save_B_vlsv(input_tuple):
    '''
    input_tuple[0]: run (string)  # 'EGL' or 'FHA'
    input_tuple[1]: fileIndex (int)
    '''
    # instantiate VlsvWriter object
    run, fileIndex = input_tuple
    filename = get_vlsvfile_fullpath( run, fileIndex)
    f = pt.vlsvfile.VlsvReader( filename )      # f contains the vg_ mesh over which Biot-Savart is integrated
    if run == 'EGL':
        f_J_sidecar = pt.vlsvfile.VlsvReader('/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/visualizations_2/ballooning/jlsidecar_bulk1.egl.{}.vlsv'.format(str(fileIndex).zfill(7)))
        f_iono = pt.vlsvfile.VlsvReader( '/wrk-vakka/group/spacephysics/vlasiator/temp/ionogrid_FHA.vlsv' )
        save_dir = '/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/sidecars/ig_B/'
    elif run == 'FHA':
        f_J_sidecar = None
        f_iono = f
        save_dir = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/ig_B/'
    elif run == 'FIA':
        f_J_sidecar = None
        f_iono = f
        save_dir = '/wrk-vakka/group/spacephysics/vlasiator/3D/FIA/bulk_sidecars/ig_B/'
    # calculate magnetic fields
    ig_r = get_ig_r(f_iono)                     # f_iono contains the ionospheric mesh (the locations where B is evaluated)
    B_iono = B_ionosphere(f, ig_r = ig_r, method = "integrate")
    try:   # UNTESTED (as of 28.7.2023 the .vslv sidecar and visit plots assumed r_io=5 RE for EGL, FHA, FIA. need to rerun with this try/except block )
        r_io = float(f.get_config()['ionosphere']['downmapRadius'][0]) * R_EARTH
    except: #EGL
        r_io = 5 * 6.371e6
    B_inner, B_outer = B_magnetosphere(f, f_J_sidecar = f_J_sidecar, r_io = r_io, ig_r = ig_r)
    # write to file
    filename_vlsv = save_dir + 'ionosphere_B_sidecar_{}.{}.vlsv'.format(run, str(fileIndex).zfill(7))
    mkdir_path(filename_vlsv)
    writer = pt.vlsvfile.VlsvWriter(f_iono, filename_vlsv)
    writer.write(ig_r,'ig_r','VARIABLE','ionosphere')
    writer.write(B_iono,'ig_B_ionosphere','VARIABLE','ionosphere')
    writer.write(B_inner,'ig_B_inner','VARIABLE','ionosphere')
    writer.write(B_outer,'ig_B_outer','VARIABLE','ionosphere')
    return ig_r, B_iono, B_inner, B_outer



def plot_Dst(run, start, stop, step):
    t = np.arange(start, stop +1, step )
    #Dsts = np.zeros(int((stop - start + 1) / step))
    Dsts = np.array(list(range(start, stop, step)))*0.0   # weird way to initialize array, but want the for loop to work
    for i, fileIndex in enumerate(range(start, stop, step)):
    #for i, fileIndex in enumerate(range(start, stop + 1, step)):
        try:
            if run == 'EGL':
                f_J_sidecar = pt.vlsvfile.VlsvReader('/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/visualizations_2/ballooning/jlsidecar_bulk1.egl.{}.vlsv'.format(str(fileIndex).zfill(7)))
            else:
                f_J_sidecar = None
            filename = get_vlsvfile_fullpath(run, fileIndex)
            f = pt.vlsvfile.VlsvReader(filename)
            Dst = calc_Dst(f, f_J_sidecar = f_J_sidecar, r_io = 5 * 6.371e6)
            Dsts[i] = Dst
        except:   # IndexError?
            print("Tried index {}, Dsts array size {}".format(i, Dsts.size))
    plt.plot(t, Dsts * 1e9)
    plt.xlabel('t [sec]')
    plt.ylabel('Dst [nT]')
    figfile = '/wrk-vakka/users/horakons/carrington/plots/{}/Dst/Dst_{}.png'.format(run, run)
    savefile = '/wrk-vakka/users/horakons/carrington/data/{}/Dst/Dst_{}.pickle'.format(run, run)
    mkdir_path(figfile)
    mkdir_path(savefile)
    plt.savefig(figfile)
    save(savefile, t = t, Dsts = Dsts)


def plot_Dsts(runs):
    # overplot the (saved) Dst data for every run in the list of strings "runs"
    s = ''
    for run in runs:
        s = s+'_'+run
        savefile = '/wrk-vakka/users/horakons/carrington/data/{}/Dst/Dst_{}.pickle'.format(run, run)
        data = restore(savefile)
        plt.plot(data['t'], data['Dsts']*1e9, label = run)
    plt.xlabel('time [sec]')
    plt.ylabel('Dst [nT]')
    plt.legend()
    figfile = '/wrk-vakka/users/horakons/carrington/plots/{}/Dst/Dst{}.png'.format(runs[0], s)
    plt.savefig(figfile)



# add .vlsv Writer to main code for any run and time step
# TODO: Write loop that calculates Biot-Savart and writes .vlsv files
# Make it faster --- assess bottlenecks
# New code that loops through .vlsvs and calculates dB/dt and GICs (and .vlsv sidecars?)
#  - Fourier? 

if __name__ == '__main__':
    run = ARGS.run
    if run == 'EGL':
        first = 621
        last = 1760
    elif run == 'FHA':
        first = 501
        last = 1612
    elif run == 'FIA':
        first = 1
        last = 865

 
    # 1. integrate Biot-Savart and save output into .vlsv files (modify biot_savart.sh to be multi-processor)
    from multiprocessing import Pool
    pool = Pool(int(ARGS.nproc))
    start = first + (int(ARGS.task) * int(ARGS.nproc))
    stop = start + int(ARGS.nproc)
    print('start:, ', start, ', stop: ', stop)
    input_list = [(run, i) for i in range(start, stop)]
    f_out = pool.map(save_B_vlsv, input_list)
    pool.close()
    pool.join()

    '''
    # 2. (requires #1) make a plot of Dst vs time (modify biot_savart.sh to use 1 node)
    plot_Dst(run, first, last, 20)   # FIA plot wasn't working due to file permissions
    
    # 3. (requires #1 and #2) make Dst plots as in #2, of multiple runs
    plot_Dsts(['EGL', 'FHA', 'FIA'])
    '''






