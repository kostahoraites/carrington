import numpy as np
import pytools as pt   # pytools is an Analysator package
import scipy.io

def fg_grid(f, fg_b = None):
    # get the face positions of the cells
    # note that e.g. xmin, xmax give the leftmost and rightmost cell BOUNDARIES (the volumes are between these boundaries)
    if (fg_b is None):
        fg_b = f.read_fsgrid_variable('fg_b')        # EGL: fg_b.shape = (1024, 736, 736, 3)
    dx = (f.read_parameter('xmax') - f.read_parameter('xmin')) / fg_b.shape[0]
    dy = (f.read_parameter('ymax') - f.read_parameter('ymin')) / fg_b.shape[1]
    dz = (f.read_parameter('zmax') - f.read_parameter('zmin')) / fg_b.shape[2]
    x = np.linspace( f.read_parameter('xmin') + dx/2, f.read_parameter('xmax') - dx/2, fg_b.shape[0] )
    y = np.linspace( f.read_parameter('ymin') + dy/2, f.read_parameter('ymax') - dy/2, fg_b.shape[1] )
    z = np.linspace( f.read_parameter('zmin') + dz/2, f.read_parameter('zmax') - dz/2, fg_b.shape[2] )
    return x, y, z


#indices = [621, 1760]
indices = [1491]
for index in indices:
    #vlsvfile = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{}.vlsv".format( str(index).zfill(7) )
    vlsvfile = "/wrk-vakka/users/ykempf/ionosphere/EGI/FAC_fgb/bulk.{}.vlsv".format( str(index).zfill(7) )
    f = pt.vlsvfile.VlsvReader(vlsvfile)
    x1d, y1d, z1d = fg_grid(f)
    xmin = 604; xmax = 804
    ymin = 268; ymax = 468
    zmin = 268; zmax = 468
    x1d = x1d[xmin:xmax]; y1d = y1d[ymin:ymax]; z1d = z1d[zmin:zmax]
    X, Y, Z = np.meshgrid(x1d, y1d, z1d, indexing='ij')
    fg_b = f.read_variable('fg_b')
    Bx = fg_b[xmin:xmax,ymin:ymax,zmin:zmax,0]
    By = fg_b[xmin:xmax,ymin:ymax,zmin:zmax,1]
    Bz = fg_b[xmin:xmax,ymin:ymax,zmin:zmax,2]
    #matfile = 'fg_b_EGL_index_{}.mat'.format(index)
    matfile = 'fg_b_EGILIKE_index_{}.mat'.format(index)
    scipy.io.savemat(matfile, mdict={'x1d':x1d, 'y1d':y1d, 'z1d':z1d, 'X':X, 'Y':Y, 'Z':Z, 'Bx':Bx, 'By':By, 'Bz':Bz})     # in Matlab, use dlmread() to read the .txt file?

# Now load in the data from the .mat that was just saved
#matdata = scipy.io.loadmat(matfile)



