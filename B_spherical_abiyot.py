import numpy as np
import pytools as pt

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


def cartesian_to_spherical_vector(vx, vy, vz, x, y, z):
    '''
    Convert cartesian vector(s) with coordinates (vx, vy, vz)
    at the position(s) theta, phi (note: position r does not affect the vector transformation)
    to spherical coordinates (v_r, v_theta, v_phi)

    dimensions of vx, vy, vz, x, y, z arrays must either match or be a single number
    '''
    r, theta, phi = cartesian_to_spherical(x, y, z)
    v_r =     vx * np.sin(theta) * np.cos(phi) + vy * np.sin(theta) * np.sin(phi) + vz * np.cos(theta)
    v_theta = vx * np.cos(theta) * np.cos(phi) + vy * np.cos(theta) * np.sin(phi) - vz * np.sin(theta)
    v_phi =  -vx * np.sin(phi) + vy * np.cos(phi)
    return v_r, v_theta, v_phi


f = pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/ig_B/ionosphere_B_sidecar_FHA.0001001.vlsv")

ig_r = f.read_variable('ig_r')
ig_B_ionosphere = f.read_variable('ig_B_ionosphere')
ig_B_inner = f.read_variable('ig_B_inner')
ig_B_outer = f.read_variable('ig_B_outer')

ig_B = ig_B_ionosphere + ig_B_inner + ig_B_outer

B_r, B_theta, B_phi = cartesian_to_spherical_vector(ig_B[:,0], ig_B[:,1], ig_B[:,2], ig_r[:,0], ig_r[:,1], ig_r[:,2])
