import numpy as np
from scipy import interpolate
from scipy.integrate import odeint
import ftest as ft
from myutils import timer, save, restore, get_vlsvfile_fullpath
from memory_profiler import profile

import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Tkagg')

global R_EARTH
R_EARTH = 6.371e6

# Input parameters
import argparse
parser = argparse.ArgumentParser()

parser.add_argument('-nproc', default=1, help="number of processors to use " )
global ARGS
ARGS = parser.parse_args()
#nproc = ARGS.nproc

global filename
global vlsvReader
global interpolators

def interpolators_fg(f):
    '''
    make E and B interpolators
    works for runs with field grid (fg_) data. e.g., runs  EGI, EGL, EGP (occasional FHA snapshots?)
    '''
    # Read cellids in order to sort variables
    cellids = f.read_variable("CellID")
    xmin = f.read_parameter('xmin')
    xmax = f.read_parameter('xmax')
    ymin = f.read_parameter('ymin')
    ymax = f.read_parameter('ymax')
    zmin = f.read_parameter('zmin')
    zmax = f.read_parameter('zmax')

    # Read face_B:
    face_B = f.read_variable('fg_b')
    # Read face_E:
    face_E = f.read_variable('fg_e')

    xsize = face_B.shape[0]
    ysize = face_B.shape[1]
    zsize = face_B.shape[2]

    sizes = np.array([xsize, ysize, zsize])
    maxs = np.array([xmax, ymax, zmax])
    mins = np.array([xmin, ymin, zmin])
    dcell = (maxs - mins)/(sizes.astype('float'))

    # Create x, y, and z coordinates:
    x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
    y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
    z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]

    # Create grid interpolation
    intp_E_x = interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), face_E[:,:,:,0])
    intp_E_y = interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), face_E[:,:,:,1])
    intp_E_z = interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), face_E[:,:,:,2])

    intp_B_x = interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), face_B[:,:,:,0])
    intp_B_y = interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), face_B[:,:,:,1])
    intp_B_z = interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), face_B[:,:,:,2])

    E = lambda x: np.array([intp_E_x(x), intp_E_y(x), intp_E_z(x)]).flatten()
    B = lambda x: np.array([intp_B_x(x), intp_B_y(x), intp_B_z(x)]).flatten()
    
    return E, B




def lorentz(X, t, q_over_m, E, B, time_sign):
        """
        (credit to Christian Hill)
        The equations of motion for the Lorentz force on a particle with
        q/m given by q_over_m. X=[x,y,z,vx,vy,vz] defines the particle's
        position and velocity at time t: F = ma = (q/m)[E + v×B].
        
        """
        x = X[0:3] 
        v = X[3:]
        try:
            drdt = time_sign * v
            dvdt = time_sign * q_over_m * (E(x) + np.cross(v, B(x)))
        except:     # if trace would leave simulation, hold it fixed in space
            drdt = x * 0.
            dvdt = v * 0.
        return np.hstack((drdt, dvdt))



@timer
def trace_particle(f, x0, v0, interpolators = None, particle = 'electron', time_sign = 1, nt = None, dt = 1e-3):
    '''
    f: vlsvReader object
    x0: 3-element array, initial position [m]
    v0: 3-element array, initial velocity [m/s]
    coord_list: list of 3-element coordinate arrays [m]
    particle: 'electron' or 'proton'
    time_sign: 1 or -1
    dt: time step
    '''

    if particle == 'electron':
        mass = 9.10938e-31
        charge = -1.60217e-19
    elif particle == 'proton':
        mass = 1.67262e-27
        charge = 1.60217e-19

    if interpolators is None:
        E, B = interpolators_fg(f)
    else:
        E, B = interpolators

    x_out = np.zeros([nt, 3])
    v_out = np.zeros([nt, 3])

    # Initial positon and velocity components.
    X0 = np.hstack((x0, v0))
    t = np.linspace(0, nt*dt, nt+1)
    # Do the numerical integration of the equation of motion.
    print('integrating motion...')
    X = odeint(lorentz, X0, t, args=(charge/mass,E,B,time_sign,))

    x_out = X[:, 0:3]
    v_out = X[:, 3:]
    return x_out, v_out, t



def sample_evdf(vlsvReader, x, v):
    '''
    given a vlsvReader object, find the value of the distribution function f
    at a given coordinate in position/velocity phase space: f(x,v)
    assumes a maxwellian f, and polytropic equation of state (adiabatic index 5/3)

    Inputs:
    x: position vector [m], 3-element array
    v: velocity vector [m/s], 3-element array
    '''
    gamma = 5. / 3.
    try:
        # Generally, it should be possible to find parameters with vlsvReader.get_config().
        # Tested: FHA works, but not EGL
        n0 = float(f.get_config()['proton_Magnetosphere']['rho'][0])
        T0 = 1.380649e-23 * float(f.get_config()['proton_Magnetosphere']['T'][0])
    except:
        n0 = 1e6              # EGI, EGL, FHA [m^-3]. 
        T0 = 1.380649e-23 * 5e5         # EGI, EGL [Joules]
    C = T0 * n0**(1. - gamma)
    n = vlsvReader.read_interpolated_variable('proton/vg_rho', x)
    T = C * n**(gamma - 1.)
    vbulk = vlsvReader.read_interpolated_variable('proton/vg_v', x)
    v_plasma = v - vbulk  # velocity vector in plasma frame
    v_mag = np.sqrt( v_plasma[0]**2 + v_plasma[1]**2 + v_plasma[2]**2 )
    m_e = 9.10938e-31
    f = n * (m_e / (2 * np.pi * T))**1.5 * np.exp(-m_e * v_mag**2 / (2 * T))
    return f


@timer
#@profile
def liouville_1pt(input_tuple):
    #vlsvReader, x0, v0, interpolators, particle, time_sign, nt, dt = input_tuple
    x0, v0, particle, time_sign, nt, dt = input_tuple
    x, v, t = trace_particle(vlsvReader, x0, v0, interpolators = interpolators, particle = particle, time_sign = time_sign, nt = nt, dt = dt)
    print('final x: ', x[-1])
    print('final v: ', v[-1])
    f = sample_evdf(vlsvReader, x[-1], v[-1])
    #if save_data:
    #    save('/wrk-vakka/users/horakons/carrington/data/particle_tracer/paths/f_liouville_test_{}_nt_{}_v{}.pickle'.format(particle, nt, v[0]), x=x, v=v, t=t)
    return f, x[-1], v[-1], x, v, t



def liouville(x, run, fileIndex, save_data = False, particle = 'electron', time_sign = 1, nt = 10000, dt = 1e-3, vmin = -4e6, vmax = 4e6, nv = 2 ):
#def liouville(vlsvReader, x, save_data = False,  interpolators = None, particle = 'electron', time_sign = 1, nt = 10000, dt = 1e-3 ):
    # note: vlsvReader and interpolators are global variables to enable Pool.map
    B = vlsvReader.read_interpolated_variable('vg_b_vol', x)
    vbulk = vlsvReader.read_interpolated_variable('proton/vg_v', x)
    #vbulk_mag = vbulk / np.sqrt( vbulk[0]**2 + vbulk[1]**2 + vbulk[2]**2 )
    vpar_hat = B / np.sqrt( B[0]**2 + B[1]**2 + B[2]**2 )
    vperp_hat = np.cross(np.array([1,0,0]), vpar_hat )
    vperp_hat = vperp_hat / np.sqrt( vperp_hat[0]**2 + vperp_hat[1]**2 + vperp_hat[2]**2 )
    dv = (vmax - vmin) / nv
    vpar, vperp = np.meshgrid( np.arange(vmin, vmax, dv), np.arange(vmin,vmax, dv), indexing='ij')
    from multiprocessing import Pool
    pool = Pool(int(ARGS.nproc))
    xv_list = []        # initial conditions
    for i in range(nv):
        for j in range(nv):
            xv_list.append( (x, vbulk + (vpar_hat * vpar[i, j]) + (vperp_hat * vperp[i, j]), particle, time_sign, nt, dt) )
            #xv_list.append( (x, (vpar_hat * (vbulk_mag + vpar[i, j])) + vperp_hat * vperp[i, j], particle, time_sign, nt, dt) )
            #xv_list.append( (x, (vpar_hat * (vbulk_mag + vpar[i, j])) + vperp_hat * vperp[i, j], interpolators, particle, time_sign, nt, dt) )
    with Pool(int(ARGS.nproc)) as p:
        data = p.map(liouville_1pt, xv_list)
    f = np.zeros(vpar.shape)
    x_f = np.zeros([vpar.shape[0], vpar.shape[1], 3])
    v_f = np.zeros([vpar.shape[0], vpar.shape[1], 3])
    x_trace = [[] for _ in range(nv)]     # DON'T use [[]] * nv (each list will have the same address)
    v_trace = [[] for _ in range(nv)]     # DON'T use [[]] * nv (each list will have the same address)
    t = [[] for _ in range(nv)]     # DON'T use [[]] * nv (each list will have the same address)
    for i in range(nv):
        for j in range(nv):
            ind = i*nv + j
            f[i, j] = data[ind][0]
            x_f[i, j,:] = data[ind][1]
            v_f[i, j,:] = data[ind][2]
            x_trace[i].append(data[ind][3])
            v_trace[i].append(data[ind][4])
            t[i].append(data[ind][5])
    if save_data:
        save('/wrk-vakka/users/horakons/carrington/data/particle_tracer/f_liouville_test_{}_{}_{}_nt_{}_x{:.1f}_y{:.1f}_z{:.1f}.pickle'.format(run, fileIndex, particle, nt, x[0]/R_EARTH, x[1]/R_EARTH, x[2]/R_EARTH), run = run, fileIndex = fileIndex, f = f, vpar = vpar, vperp = vperp, x_f = x_f, v_f = v_f, x = x_trace, v=v_trace, t=t, filename = filename)
    return f, vpar, vperp, x_f, v_f, x_trace, v_trace, t


#implement:
# save other variables: full x(t), v(t) for all vpar, vper pairs (keyword flag?)
# mach_f, T_f, vsw_f, n_f   to see if final point is in the solar wind
# dynamically decide how many iterations to run? Either at th outset, some kind of while loop, or try/except
# some way to break the distribution into multiple 8x8 blocks
# positrons?


if __name__ == '__main__':
    #run = 'EGL'
    #fileIndex = 1760
    run = 'EGI'
    fileIndex = 1506
    filename = get_vlsvfile_fullpath( run, fileIndex)
    vlsvReader = ft.f(filename)
    interpolators = interpolators_fg(vlsvReader)
    # "initial" conditions for (back-)tracing
    x0 = R_EARTH * np.array([11.0,0,0])
    particle = 'electron'
    nt = 80000
    #nt = 600
    dt = 5e-4           # estimate: omega_p ~ 1 Hz, omega_e ~ 1 kHz in solar wind
    time_sign = -1
    # trace particle trajectories
    #v0 = np.array([0,0,0])
    #x, v, t = trace_particle(f, x0, v0, interpolators = interpolators, particle = particle, time_sign = time_sign, nt = nt, dt = dt)
    if particle == 'electron': 
        vmin = -1e7
        vmax = 1e7
    elif particle == 'proton':
        vmin = -1e5
        vmax = 1e5
    nv = 8
    f, vpar, vperp, x_f, v_f, x, v, t = liouville(x0, run, fileIndex, save_data = True, particle = particle, time_sign = time_sign, nt = nt, dt = dt, vmin = vmin, vmax = vmax, nv = nv )
    print(vpar)
    print(vperp)
    print(f)
    # make plot
    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #ax.plot3D(x[:,0]/R_EARTH, x[:,1]/R_EARTH, x[:,2]/R_EARTH)
    #ax.set_xlabel(r'x [$R_E$]')
    #ax.set_ylabel(r'y [$R_E$]')
    #ax.set_zlabel(r'z [$R_E$]')
    #plt.savefig('../test/particle_tracer_nt{}_{}.png'.format(nt, particle))


