import pytools as pt
import numpy as np
from multiprocessing import Pool
from myutils import get_vlsvfile_fullpath, save

global R_EARTH
R_EARTH = 6378137.0  # 6.371e6

def lorentziator_ic(vlsvobj, x_i_list, vmin, vmax, nv, fileout = "input.txt"):
    '''
    generate initial conditions: a regular vpar,vper grid in the range [vmin, vmax]
    x_i_list: list of 3-element coordinates (array or list) to initialize particles
    
    outputs: the initial position and velocity coordinates, x_i and v_i
    fileout: text file containing the comma-separated coordinates
    '''
    x_i = np.zeros([len(x_i_list), nv, nv, 3])
    v_i = np.zeros([len(x_i_list), nv, nv, 3])
    dv = (vmax - vmin) / (nv-1)
    vpar, vperp = np.meshgrid( np.arange(vmin, vmax+dv, dv), np.arange(vmin,vmax+dv, dv), indexing='ij')
    f = open(fileout, "w")

    nx = len(x_i_list)
    for ix, x in enumerate(x_i_list):
        B = vlsvobj.read_interpolated_variable('vg_b_vol', x)
        vbulk = vlsvobj.read_interpolated_variable('proton/vg_v', x)
        vpar_hat = B / np.sqrt( B[0]**2 + B[1]**2 + B[2]**2 )
        vperp_hat = np.cross(np.array([1,0,0]), vpar_hat )
        vperp_hat = vperp_hat / np.sqrt( vperp_hat[0]**2 + vperp_hat[1]**2 + vperp_hat[2]**2 )
        for i in range(nv):
            for j in range(nv):
                v = vbulk + (vpar_hat * vpar[i, j]) + (vperp_hat * vperp[i, j])
                if (i == (nv-1)) & (j == (nv-1)) & (ix == nx):
                    suffix = ""
                else:
                    suffix = "\n"
                f.write("{},{},{},{},{},{}".format(x[0], x[1], x[2],
                                                   v[0], v[1], v[2])+suffix)
                x_i[ix, i, j, :] = x
                v_i[ix, i, j, :] = v

    f.close()
    return x_i, v_i



if __name__ == "__main__":
    vmin = -5e6     # ~1keV
    vmax = 3.75e6
    nv = 8        # phase space nv^2 elements (vpar x vper)
    print('') 
    run = 'EGI'
    fileIndex = 1300  # time to start tracing.
    vlsvobj = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(run, fileIndex))
    #x_i_list = [ [(11.5 + i*1.5)* R_EARTH, 0*R_EARTH, 0*R_EARTH] for i in range(2)]
    x_i_list = [ [13* R_EARTH, 0*R_EARTH, 0*R_EARTH] ]
    fileout = "input.txt"
    x_i, v_i = lorentziator_ic(vlsvobj, x_i_list, vmin, vmax, nv, fileout = fileout)
    # save metadata
    # note time step information needs to be input by hand separately. It's a bit klugy but just saving it here in the metadata file
    dt = -0.5
    t = np.arange(1300, 1298.5, dt)
    save("metadata.pck", vmin = vmin, vmax = vmax, nv = nv, run = run, fileIndex = fileIndex, x_i_list = x_i_list, x_i = x_i, v_i = v_i, t = t, dt = dt)
    #fileout = "/wrk-vakka/users/horakons/carrington/data/particle_tracer/input_lorentziator.txt"
    
