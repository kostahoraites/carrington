import matplotlib.pyplot as plt
import argparse
import pytools as pt
import numpy as np
import scipy
from pyPlots import plot_vdf

import pandas as pd

from myutils import get_vlsvfile_fullpath

from carrington_beta_star import fit_magnetopause

global run
run = 'FHA'   # 'EGI','EGL' 'EGP', 'FHA'

def parallelize_this(t):
    # input t is the time in seconds [integer]
    # returns: magnetopause position at time t [float]
    # FILL THIS IN
    
    fileIndex = str(t) #time in seconds

    if len(fileIndex) == 3:
        filenumber = '0000'+fileIndex
    if len(fileIndex) == 4:
        filenumber = '000'+fileIndex

    filename = get_vlsvfile_fullpath(run, fileIndex)
    f = pt.vlsvfile.VlsvReader(filename) # vlsvReaderObject , from reading the .vlsv file
    #f = pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{}.vlsv".format(filenumber)) # vlsvReaderObject , from reading the .vlsv file

    #root_dir = '/wrk-vakka/users/emilirin/plots/' # where you want to save the plot
    root_dir = '/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/new_plots/'


    threshold = 0.5  # beta_star
    delta = 0.05     # consider data satisfying beta_star = threshold +/- delta
                     # no use for delta anymore!

    #mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = False)  # magnetopause position at time t 
    #original way of saving the data
    #with open('mp_nose_x.txt', 'a') as f:
    #    f.write(str(mp_nose_x))
    #    f.write(' ')
    #    f.write(str(t))
    #    f.write('\n')

    #f=pt.vlsvfile.VlsvReader(filename)
    #f.optimize_open_file()
    #t=f.read_parameter("time")
    #if t == None:
    #   t=f.read_parameter("t")
    #if t == None:
    #       print("Unknown time format in file " + filename)

    mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = True)  # magnetopause position at time t
    RE = 6371000        # m
    #dx = 4000000       # look behind the magnetopause a short distance --- 4 cells?
    Rmax = mp_nose_x + 5.  # integrate density out to this distance (in RE) to estimate inertia term C_I (see Freeman 1998)
    dx = RE / 2         # look behind the magnetopause a short distance --- 4 cells?
    coord_mp = np.array([(mp_nose_x * RE)-dx, 0, 0])
    coord_bs = np.array([Rmax*RE, 0, 0])
    #cellid = f.get_cellid(coord)
    #B_list = f.read_variable('vg_b_vol',operator='magnitude',cellids=cellid)
    B = f.read_interpolated_variable('vg_b_vol', coord_mp, operator='magnitude')
    n = f.read_interpolated_variable('proton/vg_rho', coord_bs)
    P_d = f.read_interpolated_variable('proton/vg_pdyn', coord_bs)
    v = f.read_interpolated_variable('proton/vg_v', coord_bs, operator='magnitude')
    #Rmax = 20  # integrate density out to this distance (in RE) to estimate inertia term C_I (see Freeman 1998)
    step = 0.1
    rs = RE * np.arange(mp_nose_x - (dx/RE),Rmax,step)
    coords=list(np.transpose( np.array([rs, np.zeros(rs.size), np.zeros(rs.size)])))
    rhos = np.zeros(len(coords))
    m_proton = 1.67e-27
    for i, temp_coord in enumerate(coords):
        rhos[i] = f.read_interpolated_variable('proton/vg_rho', temp_coord) * m_proton
    C_I = np.sum(np.array(rhos)) * step * RE     # ~ \int_x1^Rmax rho dx, where x1 is magnetopause position ( integrated density across sheath )
    #return (mp_nose_x, B)
    return {'B':B, 'mp_nose_x':mp_nose_x, 'C_I':C_I, 'n':n, 'P_d':P_d, 'v':v}




#Parse command line arguments:
parser = argparse.ArgumentParser()
parser.add_argument('-nproc', default=1, help="number of processes, multithreading" )
args = parser.parse_args()


if run == 'EGL':
    times = np.array(np.arange(621, 1761))
elif run == 'EGI':
    times = np.array(np.arange(662, 1507))
elif run == 'EGP':
    times = np.array(np.arange(352, 506))    #or 269-506. But doesn't exist precipitation data before t=352
elif run == 'FHA':
    times = np.array(np.arange(501, 1498))
B = times * 0.
standoff_dist = times * 0.
C_I = times * 0.
n = times * 0.
P_d = times * 0.
v = times * 0.



## Parallel processing
from multiprocessing import Pool
pool = Pool(int(args.nproc))
#standoff_dist = pool.map(parallelize_this, times)
#tuples = pool.map(parallelize_this, times)
dicts = pool.map(parallelize_this, times)
pool.close()
pool.join()

for i, d in enumerate(dicts):
    standoff_dist[i] = d['mp_nose_x']
    B[i] = d['B']
    C_I[i] = d['C_I']
    n[i] = d['n']
    P_d[i] = d['P_d']
    v[i] = d['v']

#for i, tup in enumerate(tuples):
#    standoff_dist[i] = tup[0]
#    B[i] = tup[1]

print('TEST FLAG')
print(times)
print(standoff_dist)

#new way of saving the data
#i_sort = np.argsort(np.array(t))
#ts = np.array(t)[i_sort]
#xs = np.array(mp_nose_x)[i_sort]
#df = pd.DataFrame(data={'x':xs, 't':ts})
#df.to_csv('mp_nose_x.csv', index=False)
df = pd.DataFrame(data={'x':standoff_dist, 'B':B, 't':times, 'C_I':C_I, 'n':n, 'P_d':P_d, 'v':v})
df.to_csv('mp_nose_x_new_{}.csv'.format(run), index=False)


# now plot the results
#plt.plot(times, standoff_dist)
#plt.title(r'Magnetopause standoff distance')
#plt.xlabel(r'Time [s]')
#plt.ylabel(r'Standoff distance [$R_{E}$]')
#plt.savefig("standoff_distance_vs_time.png")

