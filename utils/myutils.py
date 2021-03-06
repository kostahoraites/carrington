import os
import numpy as np

'''
 show() creates uses display to display a plot created by tell(),
 which is saved in a default location
 the input, plt, is assumed to be a matplotlib plot

    ex. (on Turso login node, won't work on a compute note)
 import matplotlib.pyplot as plt
 from utils import *
 plt.plot([1,2,3])
 tell(plt)
 show()

Note that show (and tell?) won't work in an interactive session.

Another (more useful) option is to just type 'show' into the command line on turso

''' 

def get_bulklocation(run, fileIndex):
    # load data
    if run.upper() == 'EGI':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/"
    elif run.upper() == 'EGL':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/"
    elif run.upper() == 'EGM':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGM/bulk/"
    elif run.upper() == 'EGN':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGM/"
    elif run.upper() == 'EGO':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGO/"
    elif run.upper() == 'EGP':
        location = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGP/bulk1/"
        #location = "/wrk-vakka/group/spacephysics/vlasiator/3D/{}/bulk5/".format(run.upper())
    elif run.upper() == 'EGILIKE':
        location = "/wrk-vakka/group/spacephysics/vlasiator/temp/EGI_like/"
        #location = "/wrk-vakka/group/spacephysics/vlasiator/3D/{}/bulk5/".format(run.upper())
    return location


def get_filename(run, fileIndex):
    # load data
    if run.upper() == 'EGI':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGL':
        filename = "bulk1.{}.{}.vlsv".format(run.lower(), str(fileIndex).zfill(7) )
    elif run.upper() == 'EGM':
        filename = "bulk.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGN':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGO':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGP':
        filename = "bulk1.{}.vlsv".format(str(fileIndex).zfill(7) )
        #filename = "bulk5.{}.vlsv".format(str(fileIndex).zfill(7) )
    elif run.upper() == 'EGILIKE':   # test run
        filename = "bulk.{}.vlsv_fg".format(str(fileIndex).zfill(7) )
    return filename



def get_vlsvfile_fullpath(run, fileIndex):
    return get_bulklocation(run, fileIndex) + get_filename(run, fileIndex)


def tell(plt):
   tmp_dir = '/wrk-vakka/users/horakons/temp/'
   plt.savefig(temp_dir+'plot.png', dpi=300, bbox_inches='tight')


def show():
   tmp_dir = '/wrk-vakka/users/horakons/temp/'
   os.system('display ' + tmp_dir + 'plot.png')


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


def mkdir_path(path):
    filedir_list = path.split('/')
    filedir = path[:-len(filedir_list[-1])]
    if not(os.path.exists(filedir)):
         os.system('mkdir -p {}'.format(filedir))


def numcurl3d(inputarray, CELLSIZE_XYZ):
    # Assumes input array is of format [nx,ny,nz,3]
    # CELLSIZE_XYZ is 3-element array or list of grid spacings
    jac = numjacobian3d(inputarray, CELLSIZE_XYZ)
    # Output array is of format [nx,ny,nz,3]
    curl = np.zeros(inputarray.shape)
    curl[:,:,:,0] = jac[:,:,:,2,1]-jac[:,:,:,1,2]
    curl[:,:,:,1] = jac[:,:,:,0,2]-jac[:,:,:,2,0]
    curl[:,:,:,2] = jac[:,:,:,1,0]-jac[:,:,:,0,1]
    return curl


def numjacobian3d(inputarray, CELLSIZE_XYZ):
    # Assumes input array is of format [nx,ny,nz,3]
    # CELLSIZE_XYZ is 3-element array or list of grid spacings
    nx,ny,nz = inputarray[:,:,:,0].shape
    jac = np.zeros([nx,ny,nz,3,3])
    jac[:,:,:,0,0], jac[:,:,:,0,1], jac[:,:,:,0,2] = np.gradient(inputarray[:,:,:,0], CELLSIZE_XYZ[0])
    jac[:,:,:,1,0], jac[:,:,:,1,1], jac[:,:,:,1,2] = np.gradient(inputarray[:,:,:,1], CELLSIZE_XYZ[1])
    jac[:,:,:,2,0], jac[:,:,:,2,1], jac[:,:,:,2,2] = np.gradient(inputarray[:,:,:,2], CELLSIZE_XYZ[2])
    # Output array is of format [nx,ny,nz,3,3]
    #  :,:,component, derivativedirection
    # so for example: dAx/dx = :,:,:,0,0
    #                 DAy/dz = :,:,:,1,2
    return jac


def save(file,**kwargs):
    """
    Save the value of some data in a file.
    Usage: 
    a = 3; bvar = 'string'; test = [True, '300', 3.14]
    save('misdatos.pickle',a=a,b=bvar,test=test)
    """
    import pickle
    f=open(file,"wb")
    pickle.dump(kwargs,f,protocol=2)
    f.close


def restore(file):
    """
    Read data saved with save function.
    Usage: datos = restore('misdatos.pickle')    #datos is a dict with the saved variables
    """
    import pickle
    f=open(file,"rb")
    result = pickle.load(f)
    f.close
    return result

