import matplotlib.pyplot as plt
import argparse
import pytools as pt
import numpy as np
import scipy
from pyPlots import plot_vdf

import pandas as pd



from carrington_beta_star import fit_magnetopause

def parallelize_this(t):
    # input t is the time in seconds [integer]
    # returns: magnetopause position at time t [float]
    # FILL THIS IN
    
    run = 'EGL'
    fileIndex = str(t) #time in seconds

    if len(fileIndex) == 3:
        filenumber = '0000'+fileIndex
    if len(fileIndex) == 4:
        filenumber = '000'+fileIndex

    f = pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{}.vlsv".format(filenumber)) # vlsvReaderObject , from reading the .vlsv file
    #root_dir = '/wrk-vakka/users/emilirin/plots/' # where you want to save the plot
    root_dir = '/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/new_plots/'

    threshold = 0.5  # beta_star
    delta = 0.05     # consider data satisfying beta_star = threshold +/- delta
                     # no use for delta anymore!

    mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = False)  # magnetopause position at time t 
    #mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = False)  # magnetopause position at time t 
    #original way of saving the data
    #with open('mp_nose_x.txt', 'a') as f:
    #    f.write(str(mp_nose_x))
    #    f.write(' ')
    #    f.write(str(t))
    #    f.write('\n')
    return mp_nose_x




#Parse command line arguments:
parser = argparse.ArgumentParser()
parser.add_argument('-nproc', default=1, help="number of processes, multithreading" )
args = parser.parse_args()

times = np.arange(621, 1761)

## Parallel processing
from multiprocessing import Pool
pool = Pool(int(args.nproc))
standoff_dist = pool.map(parallelize_this, times)

#new way of saving the data
#i_sort = np.argsort(np.array(t))
#ts = np.array(t)[i_sort]
#xs = np.array(mp_nose_x)[i_sort]
#df = pd.DataFrame(data={'x':xs, 't':ts})
#df.to_csv('mp_nose_x.csv', index=False)
df = pd.DataFrame(data={'x':np.array(standoff_dist), 't':np.array(times)})
df.to_csv('mp_nose_x.csv', index=False)


# now plot the results
#plt.plot(times, standoff_dist)
#plt.title(r'Magnetopause standoff distance')
#plt.xlabel(r'Time [s]')
#plt.ylabel(r'Standoff distance [$R_{E}$]')
#plt.savefig("standoff_distance_vs_time.png")

pool.close()
pool.join()
