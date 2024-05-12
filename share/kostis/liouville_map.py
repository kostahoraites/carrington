import matplotlib
import matplotlib.pyplot as plt
import ptrReader
import pytools as pt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from copy import deepcopy
import numpy as np

R_EARTH = 6.371e6

def sample_evdf(v_f):
    '''
    let initial position and velocity of a particle be x_i, v_i
    and final position and velocity be x_f, v_f

    sample_evdf(v_f) finds the value of the distribution function f(x_i, v_i)

    Liouville's theorem says distribution is constant along trajectories: 
        f(x_i, v_i) = f(x_f, v_f)
    Here, Liouville mapping assumes a (drifting) maxwellian f, 
    and polytropic equation of state (adiabatic index 5/3)

    Inputs:
    v_f: 'final' velocity vector [m/s], 3-element array
    '''
    # EGI parameters:
    gamma = 5. / 3.
    n0 = 1e6              
    T0 = 1.380649e-23 * 5e5         # [Joules]
    vbulk = np.array([-7.5e5, 0, 0])
    C = T0 * n0**(1. - gamma)
    v_plasma = v_f - vbulk  # velocity vector in plasma frame
    v_mag = np.sqrt( v_plasma[0]**2 + v_plasma[1]**2 + v_plasma[2]**2 )
    m_e = 9.10938e-31
    f = n0 * (m_e / (2 * np.pi * T0))**1.5 * np.exp(-m_e * v_mag**2 / (2 * T0))
    return f



def liouville_map(ptr_i, ptr_f):
    '''
    Inputs:
     ptr_i ("initial" .ptr file, in the magnetosheath)
     ptr_f ("final" .ptr file, in the solar wind)
    '''
    import ptrReader
    p_i = ptrReader.read_ptr_file(ptr_i)
    p_f = ptrReader.read_ptr_file(ptr_f)
    nv = int( (p_i.pos().size /3.)**0.5 )
    x_i = p_i.pos().reshape([nv,nv,3])   # shape (nv^2, 3) -> (nv,nv,3)
    v_i = p_i.vel().reshape([nv,nv,3])   # shape (nv^2, 3) -> (nv,nv,3)
    x_f = p_f.pos().reshape([nv,nv,3])   # shape (nv^2, 3) -> (nv,nv,3)
    v_f = p_f.vel().reshape([nv,nv,3])   # shape (nv^2, 3) -> (nv,nv,3)    

    # HARD CODED LORENTZIATOR CONFIG PARAMETERS:
    vmin = -1e7  # m/s
    vmax = 1e7

    # vpar-vper grid for plotting (note: .ptr initial velocities are arranged in even vpar-vper grid)
    dv = (vmax - vmin)/(nv-1)
    vpar, vperp = np.meshgrid( np.arange(vmin, vmax+dv, dv), np.arange(vmin,vmax+dv, dv), indexing='ij')

    # Liouville map the final positions
    f = np.zeros(vpar.shape)
    for i in range(nv):
        for j in range(nv):
            f[i,j] = sample_evdf(v_f[i, j, :])

    # Mask out magnetosheath data (simple test)
    r_f = np.linalg.norm(x_f, axis=-1)
    r_f_mask = (r_f < 17*R_EARTH)
    f_masked = deepcopy(f)
    f_masked[r_f_mask] = -1e-15  # dummy value (<0)

    fig, ax = plt.subplots(1,2)

    im = ax[0].pcolormesh(vpar, vperp, f_masked,
                   cmap='bwr', shading='auto', vmin=-3e-15, vmax=3e-15)
    ax[0].set_xlabel(r'$v_\parallel$')
    ax[0].set_ylabel(r'$v_\perp$')
    ax[0].set_title(r'$f_(v_\parallel, v_\perp)$')
    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    mycbar = fig.colorbar(im, cax=cax, orientation='vertical')

    ax[1].plot(vpar[:,nv//2], f_masked[:,nv//2], label = 'vpar')
    ax[1].plot(vperp[nv//2,:], f_masked[nv//2,:], label = 'vperp')
    ax[1].legend()
    ax[0].set_aspect('equal')
    plt.tight_layout()
    plt.savefig('liouville_map.png')


if __name__ == '__main__':
    ptr_i="/wrk-vakka/users/kpapadak/lorentziator/test_configs/population.0000091.ptr"
    ptr_f="/wrk-vakka/users/kpapadak/lorentziator/test_configs/population.0000092.ptr"
    liouville_map(ptr_i, ptr_f)
