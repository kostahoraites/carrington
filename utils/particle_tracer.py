import numpy as np
from scipy import interpolate
from scipy.integrate import odeint
import ftest as ft
from myutils import timer, save, restore, get_vlsvfile_fullpath  #, sidecar
#from memory_profiler import profile

import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Tkagg')

# for debugging:
import psutil

global R_EARTH
R_EARTH = 6378137.0   # 6.371e6

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
global lightspeed

lightspeed = 2.99792458e8 

def print_initial_coords(filename):
    # for generating input from using particle_tracer.py output .pickle files
    a = restore(filename)
    for i in range(8):
        for j in range(8):
            text = str(a['x_i'][i, j, :]).strip('[').strip(']').strip(',') +' ' + str(a['v_i'][i, j, :]).strip('[').strip(']').strip(',')
            print(' '.join(text.split()))


#@sidecar
def read_fg(vlsv, xrange = None, yrange = None, zrange = None):

    print('read_fg vlsv ', vlsv)

    # Read cellids in order to sort variables
    cellids = vlsv.read_variable("CellID")
    xmin = vlsv.read_parameter('xmin')
    xmax = vlsv.read_parameter('xmax')
    ymin = vlsv.read_parameter('ymin')
    ymax = vlsv.read_parameter('ymax')
    zmin = vlsv.read_parameter('zmin')
    zmax = vlsv.read_parameter('zmax')

    print("memory usage before reading fields:")
    process = psutil.Process()
    print(process.memory_info().rss / 1e9)  # in GB

    # Read face_B:
    face_B = vlsv.read_variable('fg_b')
    # Read face_E:
    face_E = vlsv.read_variable('fg_e')

    print("memory usage after reading fields:")
    process = psutil.Process()
    print(process.memory_info().rss / 1e9)  # in GB

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

    # define array indices, outlining a box-shaped subset of fg grid
    if xrange == None:
        ixmin = 0
        ixmax = x.size
        print('x OOPS!')
    else:
        ixmin = int((xrange[0] - xmin) // dcell[0])
        ixmax = int((xrange[1] - xmin) // dcell[0])
        print('x inds: ', ixmin, ixmax)

    if yrange == None:
        iymin = 0
        iymax = y.size
    else:
        iymin = int((yrange[0] - ymin) // dcell[1])
        iymax = int((yrange[1] - ymin) // dcell[1])
        print('y inds: ', iymin, iymax)

    if zrange == None:
        izmin = 0
        izmax = z.size
    else:
        izmin = int((zrange[0] - zmin) // dcell[2])
        izmax = int((zrange[1] - zmin) // dcell[2])
        print('z inds: ', izmin, izmax)

    # KLUG: need to add 0 to make face_E available for deletion
    # see    https://stackoverflow.com/questions/38324603/how-to-keep-a-slice-of-a-numpy-array-and-clear-the-rest-from-memory 
    face_E_indexed = face_E[ixmin:ixmax,iymin:iymax,izmin:izmax,:] + 0
    face_B_indexed = face_B[ixmin:ixmax,iymin:iymax,izmin:izmax,:] + 0

    del face_E
    del face_B

    return face_E_indexed, face_B_indexed, x, y, z, ixmin, ixmax, iymin, iymax, izmin, izmax, dcell


def interpolators_fg(f, xrange = None, yrange = None, zrange = None):
    '''
    make E and B interpolators
    works for runs with field grid (fg_) data. e.g., runs  EGI, EGL, EGP (occasional FHA snapshots?)
    xrange, yrange, zrange: set to two-element list/array of floats [min, max] to consider only a spatial subset of the fg grid
        e.g. xrange = [-1e8, 1e8]  # meters
    '''

    face_E_indexed, face_B_indexed, x, y, z, ixmin, ixmax, iymin, iymax, izmin, izmax, dcell = read_fg(f, xrange = xrange, yrange = yrange, zrange = zrange)

    print("memory usage after face_E,B deletion:")
    process = psutil.Process()
    print(process.memory_info().rss / 1e9)  # in GB

    # edge-centered E-field interpolation
    intp_E_x = interpolate.RegularGridInterpolator((x[ixmin:ixmax], 
                                                    y[iymin:iymax]-0.5*dcell[1], 
                                                    z[izmin:izmax]-0.5*dcell[2]), 
                                                    face_E_indexed[:,:,:,0])
    intp_E_y = interpolate.RegularGridInterpolator((x[ixmin:ixmax]-0.5*dcell[0], 
                                                    y[iymin:iymax], 
                                                    z[izmin:izmax]-0.5*dcell[2]), 
                                                    face_E_indexed[:,:,:,1])
    intp_E_z = interpolate.RegularGridInterpolator((x[ixmin:ixmax]-0.5*dcell[0],
                                                    y[iymin:iymax]-0.5*dcell[1],
                                                    z[izmin:izmax]), 
                                                    face_E_indexed[:,:,:,2])

    # face-centered B-field interpolation
    intp_B_x = interpolate.RegularGridInterpolator((x[ixmin:ixmax]-0.5*dcell[0],
                                                    y[iymin:iymax],
                                                    z[izmin:izmax]),
                                                    face_B_indexed[:,:,:,0])
    intp_B_y = interpolate.RegularGridInterpolator((x[ixmin:ixmax],
                                                    y[iymin:iymax]-0.5*dcell[1],
                                                    z[izmin:izmax]),
                                                    face_B_indexed[:,:,:,1])
    intp_B_z = interpolate.RegularGridInterpolator((x[ixmin:ixmax],
                                                    y[iymin:iymax],
                                                    z[izmin:izmax]-0.5*dcell[2]),
                                                    face_B_indexed[:,:,:,2])

    print("memory usage after making interpolators:")
    process = psutil.Process()
    print(process.memory_info().rss / 1e9)  # in GB

    E = lambda x: np.array([intp_E_x(x), intp_E_y(x), intp_E_z(x)]).flatten()
    B = lambda x: np.array([intp_B_x(x), intp_B_y(x), intp_B_z(x)]).flatten()

    print("memory usage after E, B lambda creation:")
    process = psutil.Process()
    print(process.memory_info().rss / 1e9)  # in GB
    
    return E, B


def intp_tdepify(intp_list, t_list):
    '''
     take two spatial interpolators at times t0 and t1
     and construct a time-dependent interpolator which allows time to be specified with keyword 't'
    '''
    # sort the inputs by time
    i_sort = np.argsort(np.array(t_list))
    intp_list = [intp_list[i] for i in i_sort]
    t_list = [t_list[i] for i in i_sort]

    if len(list(t_list)) != len(list(intp_list)):
        print("input lists of times and interpolators should have equal length!")
        return

    if len(t_list) == 1:
        # wrapper for the interpolator that includes dummy keyword t
        def dummy_func(x, t=None):
            return intp_list[0](x)  # keyword t doesn't do anything. (no time interpolation)
        return dummy_func
    elif len(t_list) == 2:
        # interpolate across 2 times
        intp0 = intp_list[0]
        intp1 = intp_list[1]
        t0 = t_list[0]
        t1 = t_list[1]
        print('making a two-time interpolator: ', t0, ' to ', t1)
        if t0==t1:
            print('ERROR! inputs t0 and t1 cannot be equal')
            return None
        else: #t1 > t0
            def t_func(x, t=0):
                if t < t0 or t > t1:
                    print('linear interpolation: time t out of bounds!')
                    print('t=', t, 't0=', t0, 't1=', t1)
                    return None
                else:
                    # linear time interpolation:
                    return intp1(x) * (t-t0)/(t1-t0) + intp0(x) * (t1-t)/(t1-t0)
            return t_func
        '''
        elif t0 > t1:
            # recursion: reorder the times so t1 > t0
            return intp_tdepify([intp1, intp0], [t1, t0]) 
        '''
    else: # len(t_list) > 2 
        print('tdepify with a list!')
        # construct a list of time-interpolators
        n_intervals = len(t_list)-1  # number of time intervals
        intps_t = []
        for i in range(n_intervals):
            # recursion: input lists with length 2 (see above elif case)
            intps_t.append( intp_tdepify(intp_list[i:i+2], t_list[i:i+2]) )
        t0 = t_list[0]
        delta_t = t_list[1] - t_list[0]  # approximate time difference throughout t_list
        print('t0 = ', t0)
        print('t_f = ', t_list[-1])

        # wrapper function that links together many time-interpolators
        def output_func(x, t = 0):
            #index = min( int( np.abs((t - t0) // delta_t) ), n_intervals-1 )  # initial guess (should be ~1 off at worst). klug?? this np.abs should not be necessary
            index = int( np.abs((t - t0) // delta_t) ) % n_intervals  # initial guess (should be ~1 off at worst for near-uniformly spaced times). klug?? this np.abs should not be necessary
            tdiff_sign = int(np.sign(t - t_list[index]))
            count = 0
            if tdiff_sign != 0:
                try:
                    while ((t<t_list[index]) | (t>t_list[index+1])) & (count < n_intervals):
                        #print('old index {}'.format(index))
                        # iterate (if necessary) until correct time range is found
                        #index += tdiff_sign
                        index = (index + tdiff_sign) % n_intervals
                        #print('new index {}'.format(index))
                        #print('new times t0 = {}, t1 = {}'.format(t_list[index], t_list[index+1]))
                        count += 1
                    if count == n_intervals:
                        print("output_func: count exceeded n_intervals!")
                        print("index = {}, t = {}, tdiff_sign = {}, t0 = {}, t_f = {}, delta_t = {}, n_intervals = {}, count = {}".format(index, t, tdiff_sign, t0, t_list[-1], delta_t, n_intervals, count))
                except IndexError:
                    #print("output_func: time t out of bounds!")
                    print("index = {}, t = {}, tdiff_sign = {}, t0 = {}, t_f = {}, delta_t = {}, n_intervals = {}, count = {}".format(index, t, tdiff_sign, t0, t_list[-1], delta_t, n_intervals, count))
                    import sys
                    sys.exit()
            intp_t = intps_t[index]
            return intp_t(x, t=t)
    
        return output_func



def interpolators_fg_tdep(f_list, xrange = None, yrange = None, zrange = None):
    '''
    f_list: a 2-element list of vlsvReader objects (e.g. at neighboring times)

    returns: E, B space+time interpolators (tuple)
        output interpolators are similar to the output of interpolators_fg(), defined above,
        except they also interpolate to the time t (accepted as a keyword)
    '''
    L = len(f_list)
    if L == 2:
        f0 = f_list[0]
        f1 = f_list[1]
        t0 = f0.read_parameter('time')
        t1 = f1.read_parameter('time')
        E0, B0 = interpolators_fg(f0, xrange = xrange, yrange = yrange, zrange = zrange)
        E1, B1 = interpolators_fg(f1, xrange = xrange, yrange = yrange, zrange = zrange)
        return intp_tdepify([E0, E1], [t0, t1]), intp_tdepify([B0, B1], [t0, t1])   # E and B interpolators, that also accept keyword t
    elif L > 2:
        E_list = []
        B_list = []
        t_list = []
        for i in range(L):
            print("interpolator i = {}/{}".format(i, L))
            print("memory usage:")
            process = psutil.Process()
            print(process.memory_info().rss / 1e9)  # in GB
            f0 = f_list[i]
            t0 = f0.read_parameter('time')
            E0, B0 = interpolators_fg(f0, xrange = xrange, yrange = yrange, zrange = zrange)
            E_list.append(E0)
            B_list.append(B0)
            t_list.append(t0)
        return intp_tdepify(E_list, t_list), intp_tdepify(B_list, t_list)



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
        drdt = v
        dvdt = q_over_m * (E(x, t=t) + np.cross(v, B(x, t=t)))    # time_sign deprecated (because t --> t_phys)?
        '''
        drdt = time_sign * v
        dvdt = time_sign * q_over_m * (E(x, t=t) + np.cross(v, B(x, t=t)))
        '''
    except ValueError:     # if trace would leave simulation, hold it fixed in space. See interpolate.RegularGridInterpolator(bounds_error=True) documentation
        drdt = x * 0.
        dvdt = v * 0.
    return np.hstack((drdt, dvdt))


# Boris pushers, from https://stanczakdominik.github.io/posts/on-the-recent-on-the-boris-solver-in-particle-in-cell-simulations-paper/  
# helper functions
def epsilon(q_m, electric_field, timestep):
    return q_m * timestep / 2 * electric_field 

def gamma_from_velocity(velocity):
    return np.sqrt(1 - np.sum((velocity / lightspeed)**2))

def gamma_from_u(u):
    return np.sqrt(1+np.sum((u/lightspeed)**2))

def BorisA(position, u_t_minus_half, q_m, electric_field, magnetic_field, timestep):
    # Equations 3, 6, 7a, 8, 9, 5
    uminus = u_t_minus_half + epsilon(q_m, electric_field, timestep)  # Eq. 3
    magfield_norm = np.linalg.norm(magnetic_field)
    theta = q_m * timestep / gamma_from_u(uminus) * magfield_norm  # Eq. 6
        
    b = magnetic_field / magfield_norm
    
    t = np.tan(theta/2) * b # Eq. 7a
    
    uprime = uminus + np.cross(uminus, t)  # Eq. 8
    uplus = uminus + 2/(1+(t**2).sum()) * np.cross(uprime, t)  # Eq. 9
    u_t_plus_half = uplus + epsilon(q_m, electric_field, timestep) # Eq. 5
    new_position = u_t_plus_half / gamma_from_u(u_t_plus_half) * timestep + position # Eq. 1
    return new_position, u_t_plus_half 

def BorisB(position, u_t_minus_half, q_m, electric_field, magnetic_field, timestep):
    # 3, 7b, 8, 9, 5
    uminus = u_t_minus_half + epsilon(q_m, electric_field, timestep)  # Eq. 3
    
    # Eq. 7a
    t = q_m * timestep / (2 * gamma_from_u(uminus)) * magnetic_field
    
    uprime = uminus + np.cross(uminus, t)  # Eq. 8
    uplus = uminus + 2/(1+(t**2).sum()) * np.cross(uprime, t)  # Eq. 9
    u_t_plus_half = uplus + epsilon(q_m, electric_field, timestep) # Eq. 5
    new_position = u_t_plus_half / gamma_from_u(u_t_plus_half) * timestep + position # Eq. 1
    return new_position, u_t_plus_half 
    
def BorisC(position, u_t_minus_half, q_m, electric_field, magnetic_field, timestep):
    # 3, 6, 11, 12, 5
    uminus = u_t_minus_half + epsilon(q_m, electric_field, timestep)  # Eq. 3
    magfield_norm = np.linalg.norm(magnetic_field)
    theta = q_m * timestep / gamma_from_u(uminus) * magfield_norm  # Eq. 6
    
    b = magnetic_field / magfield_norm
    
    u_parallel_minus = np.dot(uminus, b) * b # Eq. 11
    uplus = u_parallel_minus + (uminus - u_parallel_minus) * np.cos(theta) + np.cross(uminus, b) * np.sin(theta) # Eq. 12
    u_t_plus_half = uplus + epsilon(q_m, electric_field, timestep) # Eq. 5
    new_position = u_t_plus_half / gamma_from_u(u_t_plus_half) * timestep + position # Eq. 1
    return new_position, u_t_plus_half 

'''  Kostis's adaptive Boris:

       //1st order;
        earth_dipole(&r1[0],&b1[0]);
        add_ext_fields_reverse(&r1[0],&b1[0],&e1[0],t,timeline,nfiles,mmapped_e,mmapped_b, field_size_in_elements);
        Boris(&r1[0], &v1[0], &e1[0], &b1[0],dt);

        //2nd order;
        Boris(&r2[0], &v2[0], &e1[0], &b1[0],dt/2.0);
        earth_dipole(&r2[0],&b2[0]);
        add_ext_fields_reverse(&r2[0],&b2[0],&e2[0],t+dt/2.0,timeline,nfiles,mmapped_e,mmapped_b, field_size_in_elements);
        Boris(&r2[0], &v2[0], &e2[0], &b2[0],dt/2.0);

        double error[6]={
            100.0*fabs( (r2[0]-r1[0])/r1[0] ),
            100.0*fabs( (r2[1]-r1[1])/r1[1] ),
            100.0*fabs( (r2[2]-r1[2])/r1[2] ),
            100.0*fabs( (v2[0]-v1[0])/v1[0] ),
            100.0*fabs( (v2[1]-v1[1])/v1[1] ),
            100.0*fabs( (v2[2]-v1[2])/v1[2] )
        };
        double max_err= -1.0;
        for (int i=0; i<6;++i){
            if (error[i]>max_err){
               max_err=error[i];
            }
        }
        double new_dt=0.9 * fabs(dt) * min(max(sqrt(TOL / (2.0 * max_err)), 0.3), 2.0);
        if (max_err<TOL){
            x[tid]=r1[0];
            y[tid]=r1[1];
            z[tid]=r1[2];
            vx[tid]=v1[0];
            vy[tid]=v1[1];
            vz[tid]=v1[2];
            t+=dt;
            double rho = mag(r1[0],r1[1],r1[2]);
            if (rho<1.2*RE){
                alive[tid]=1;
            }
            if (rho>OUTTER_LIM){
                alive[tid]=2;
            }
        }
        pdt[tid]=new_dt;
'''

def integrate_boris(method, X0, t, q_m, E, B, time_sign, gyro_subdivide = 50.):    #args=(charge/mass,E,B,time_sign)):    # , rtol = 0.01*tol_def, atol=0.01*tol_def):
    x = X0[0:3]
    v = X0[3:]
    X = np.zeros([t.size, 6])
    X[0,:] = X0

    alerted = False
    for i in range(t.size-1):
        try:
            delta_t = abs(t[i+1] - t[i])
            # gyro_subdivide: Boris solver timestep divides delta_t into an integer number of smaller steps,
            # ensuring time_step <= (gyroperiod / gyro_subdivide)
            gyroperiod = abs( 2. * np.pi / (q_m * np.linalg.norm(B(x, t=t[i])) ))
            nt = int(np.ceil(gyro_subdivide * abs(delta_t) / gyroperiod))   # note: nt >= 1
            timestep = time_sign * delta_t / nt
            #print('t_start = {}'.format(t[i]))
            #print('t_out = {}'.format(t[i+1]))
            #print('timestep = {}\n\n'.format(timestep))
            for j in range(nt):
                t_eval = t[i] + j*timestep
                E_inpt = E(x, t=t_eval)
                B_inpt = B(x, t=t_eval)
                x, v = method(x, v, q_m, E_inpt, B_inpt, timestep)
                #print(' t_eval = {}'.format(t_eval))
                #print(' E_inpt = {}'.format(E_inpt))
                #print(' B_inpt = {}'.format(B_inpt))
        except ValueError: # Stop updating particles that exit domain. see interpolate.RegularGridInterpolator(bounds_error=True) documentation
            if not alerted:
                print('Particle escaped spatial domain! x = {}, v = {}'.format(x, v))
                alerted = True
        # save particle positions at the specified times
        X[i+1,0:3] = x
        X[i+1,3:] = v

    return X


@timer
def trace_particle(f, x0, v0, interpolators = None, particle = 'electron', time_sign = 1, nt = None, dt = 1e-3, method = 'odeint', intp_time = False, res = None):
    '''
    f: vlsvReader object
    x0: 3-element array, initial position [m]
    v0: 3-element array, initial velocity [m/s]
    coord_list: list of 3-element coordinate arrays [m]
    particle: 'electron' or 'proton'
    time_sign: 1 or -1
    dt: time step
    method: 'odeint', 'BorisA', 'BorisB', 'BorisC'
    intp_time: Boolean, set to interpolate wrt time
    res: numerical value, representing gyrosubdivide if method='Boris[A-C]', and [a,r]tol if method = 'odeint'
    '''

    if particle == 'electron':
        mass = 9.1093837e-31
        charge = -1.60217e-19
    elif particle == 'proton':
        mass = 1.67262e-27
        charge = 1.60217e-19

    if interpolators is None:
        if intp_time is True:
            print("Haven't implemented time-dependent interpolation inside this function yet. Use externally-defined interpolators.")
            return
        else:
            E, B = interpolators_fg(f)
    else:
        E, B = interpolators

    x_out = np.zeros([nt, 3])
    v_out = np.zeros([nt, 3])

    # Initial positon and velocity components.
    X0 = np.hstack((x0, v0))
    t0 = vlsvReader.read_parameter('time')
    t = np.linspace(0, nt*dt, num = nt+1)
    t_phys = t0 + t * time_sign
    print('t_phys[0:10]', t_phys[0:10])
    print('t_phys[-1:-10:-1]', t_phys[-1:-10:-1])
    # Do the numerical integration of the equation of motion.
    print('integrating motion...')
    print('time', t)
    print('method = {}'.format(method))

    boris_methods = {'BorisA':BorisA, 'BorisB':BorisB, 'BorisC':BorisC}
    if method == 'odeint':
        if res == None:
            res = 1.49012e-10 # 0.01 * Default tolerance. See scipy.odeint() documentation. ewt = rtol * abs(y) + atol.
        X = odeint(lorentz, X0, t_phys, args=(charge/mass,E,B,time_sign,), rtol = res, atol = res)
    else:
        X = integrate_boris(boris_methods[method], X0, t_phys, charge/mass, E, B, time_sign, gyro_subdivide = res)

    x_out = X[:, 0:3]
    v_out = X[:, 3:]
    return x_out, v_out, t


def sample_evdf(vlsvReader, x, v):
    '''
    given a vlsvReader object, find the value of the distribution function f
    at a given coordinate in position/velocity phase space: f(x,v)
    assumes a maxwellian f, and polytropic equation of state (polytropic index 5/3)

    Inputs:
    x: position vector [m], 3-element array
    v: velocity vector [m/s], 3-element array
    '''
    gamma = 5. / 3.
    #try:
    #    # Generally, it should be possible to find parameters with vlsvReader.get_config().
    #    # Tested: FHA works, but not EGL
    #    n0 = float(f.get_config()['proton_Magnetosphere']['rho'][0])
    #    T0 = 1.380649e-23 * float(f.get_config()['proton_Magnetosphere']['T'][0])
    #except:
    #    n0 = 1e6              # EGI, EGL, FHA [m^-3]. 
    #    T0 = 1.380649e-23 * 5e5         # EGI, EGL [Joules]
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
    x0, v0, particle, time_sign, nt, dt, method, intp_time, res = input_tuple
    x, v, t = trace_particle(vlsvReader, x0, v0, interpolators = interpolators, particle = particle, time_sign = time_sign, nt = nt, dt = dt, method = method, intp_time = intp_time, res = res)
    print('t[0:200]:', t[0:200])
    print('final x: ', x[-1])
    print('final v: ', v[-1])
    f = sample_evdf(vlsvReader, x[-1], v[-1])
    #if save_data:
    #    save('/wrk-vakka/users/horakons/carrington/data/particle_tracer/paths/f_liouville_test_{}_nt_{}_v{}.pickle'.format(particle, nt, v[0]), x=x, v=v, t=t)
    return f, x[-1], v[-1], x, v, t

 
def liouville(x, run, fileIndex, save_data = False, particle = 'electron', time_sign = 1, nt = 10000, dt = 1e-3, vmin = -4e6, vmax = 4e6, nv = 2, method = 'odeint', intp_time = False, res = None):
    # note: vlsvReader and interpolators are global variables to enable Pool.map
    B = vlsvReader.read_interpolated_variable('vg_b_vol', x)
    vbulk = vlsvReader.read_interpolated_variable('proton/vg_v', x)
    #vbulk_mag = vbulk / np.sqrt( vbulk[0]**2 + vbulk[1]**2 + vbulk[2]**2 )
    vpar_hat = B / np.sqrt( B[0]**2 + B[1]**2 + B[2]**2 )
    vperp_hat = np.cross(np.array([1,0,0]), vpar_hat )
    vperp_hat = vperp_hat / np.sqrt( vperp_hat[0]**2 + vperp_hat[1]**2 + vperp_hat[2]**2 )
    dv = (vmax - vmin) / nv
    vpar, vperp = np.meshgrid( np.arange(vmin, vmax, dv), np.arange(vmin,vmax, dv), indexing='ij')

    # MULTI-THREADING:
    from multiprocessing import Pool
    pool = Pool(int(ARGS.nproc))
    xv_list = []        # initial conditions
    for i in range(nv):
        for j in range(nv):
            xv_list.append( (x, vbulk + (vpar_hat * vpar[i, j]) + (vperp_hat * vperp[i, j]), particle, time_sign, nt, dt, method, intp_time, res) )
            #xv_list.append( (x, (vpar_hat * (vbulk_mag + vpar[i, j])) + vperp_hat * vperp[i, j], particle, time_sign, nt, dt) )
            #xv_list.append( (x, (vpar_hat * (vbulk_mag + vpar[i, j])) + vperp_hat * vperp[i, j], interpolators, particle, time_sign, nt, dt) )
    with Pool(int(ARGS.nproc)) as p:
        data = p.map(liouville_1pt, xv_list)
    f = np.zeros(vpar.shape)
    x_i = np.zeros([vpar.shape[0], vpar.shape[1], 3])
    v_i = np.zeros([vpar.shape[0], vpar.shape[1], 3])
    x_f = np.zeros([vpar.shape[0], vpar.shape[1], 3])
    v_f = np.zeros([vpar.shape[0], vpar.shape[1], 3])
    x_trace = [[] for _ in range(nv)]     # DON'T use [[]] * nv (each list will have the same address)
    v_trace = [[] for _ in range(nv)]     # DON'T use [[]] * nv (each list will have the same address)
    t = [[] for _ in range(nv)]     # DON'T use [[]] * nv (each list will have the same address)
    for i in range(nv):
        for j in range(nv):
            ind = i*nv + j
            f[i, j] = data[ind][0]
            x_i[i, j, :] = xv_list[ind][0]
            v_i[i, j, :] = xv_list[ind][1]
            x_f[i, j,:] = data[ind][1]
            v_f[i, j,:] = data[ind][2]
            x_trace[i].append(data[ind][3])
            v_trace[i].append(data[ind][4])
            t[i].append(data[ind][5])
    t0 = vlsvReader.read_parameter('time')
    t_phys = t0 + t * time_sign
    if res is None:
        res_string = 'dummy-res-string'
    else:
        res_string = '_res{:.2e}'.format(res)

    if save_data:
        save('/wrk-vakka/users/horakons/carrington/data/particle_tracer/f_liouville_test_{}_{}_{}_nt_{}_x{:.1f}_y{:.1f}_z{:.1f}_intptime{}_{}{}.pickle'.format(run, fileIndex, particle, nt, x[0]/R_EARTH, x[1]/R_EARTH, x[2]/R_EARTH, int(intp_time), method, res_string), run = run, fileIndex = fileIndex, f = f, vpar = vpar, vperp = vperp, x_i = x_i, v_i = v_i, x_f = x_f, v_f = v_f, x = x_trace, v=v_trace, t=t, t_phys = t_phys, time_sign = time_sign, method = method, res = res, filename = filename)
    return f, vpar, vperp, x_i, v_i, x_f, v_f, x_trace, v_trace, t


#implement:
# save other variables: full x(t), v(t) for all vpar, vper pairs (keyword flag?)
# mach_f, T_f, vsw_f, n_f   to see if final point is in the solar wind
# dynamically decide how many iterations to run? Either at the outset, some kind of while loop, or try/except
# some way to break the distribution into multiple 8x8 blocks
# positrons?


if __name__ == '__main__':

    # TEST (constant fields)
    #run = 'TEST'
    #fileIndex = 0
    #filename = '/wrk-vakka/users/horakons/carrington/data/particle_tracer/vlsv/E_0_-1_0_e-3_B_0_0_1_e-8.vlsv'

    # EGL
    #run = 'EGL'
    #fileIndex = 1760

    # EGI
    run = 'EGI'
    fileIndex = 1300  # 1199

    filename = get_vlsvfile_fullpath( run, fileIndex)

    vlsvReader = ft.f(filename)
    # "initial" conditions for (back-)tracing
    x0 = R_EARTH * np.array([12.8,0,0])          # [11.5, 0, 0]. Note to compare with Lorentziator, must be within box |x|,|y|,|z|< 20 RE
    xrange = [-30.*R_EARTH, 30.*R_EARTH]   # if None, interpolate across whole grid
    yrange = [-30.*R_EARTH, 30.*R_EARTH]
    zrange = [-30.*R_EARTH, 30.*R_EARTH]
    particle = 'electron'
    nt = 160000
    #nt = 600
    dt = 5e-4             # 0.1?,   estimate: omega_p ~ 1 Hz, omega_e ~ 1 kHz in solar wind
    time_sign = -1
    # trace particle trajectories
    #v0 = np.array([0,0,0])
    #x, v, t = trace_particle(f, x0, v0, interpolators = interpolators, particle = particle, time_sign = time_sign, nt = nt, dt = dt)
    if particle == 'electron': 
        vmin = -5e6 #-2e7
        vmax = 5e6 #2e7
    elif particle == 'proton':
        vmin = -1e5
        vmax = 1e5
    nv = 8

    intp_time = True
    #interpolators = interpolators_fg(vlsvReader)
    if intp_time:
        n_sec = int(nt * dt) + 1
        vlsvReader_list = [ ft.f(get_vlsvfile_fullpath(run, fileIndex + time_sign * i)) for i in range(-1, n_sec+2) ]     # input the vlsv file
        interpolators = interpolators_fg_tdep(vlsvReader_list, xrange = xrange, yrange =yrange, zrange = zrange )
        print('(time-dependent) Interpolators constructed!')
    else:
        interpolators = interpolators_fg(vlsvReader)
        interpolators = (intp_tdepify([interpolators[0]], [fileIndex]), intp_tdepify([interpolators[1]], [fileIndex]))

    method = 'BorisA'  #  'odeint', 'BorisA', 'BorisB', 'BorisC'
    print('method: ', method)
    #res = 16    # numerical value, representing gyrosubdivide in integrate_boris() or [a-r]tol in odeint()
    for i in range(9, 10):
        res = 2**(i+1)         # res = gyro_subdivide
        if method == 'odeint':
            res = res * 1e-10  # res = tol
        f, vpar, vperp, x_i, v_i, x_f, v_f, x, v, t = liouville(x0, run, fileIndex, save_data = True, particle = particle, time_sign = time_sign,
                                                            nt = nt, dt = dt, vmin = vmin, vmax = vmax, nv = nv, method = method, intp_time = intp_time, res = res)
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


