import pytools as pt
import matplotlib.pyplot as plt
import numpy as np
import scipy
from pyPlots import plot_vdf
#------------------------------------------------------------
from carrington_beta_star import fit_magnetopause
from myutils import get_vlsvfile_fullpath

#Parameters used in Horaites et al. 2023
#run = 'EGL'
#fileIndex = '857' #time in seconds

#Parameters for 2023 RCF postdoc fellow application
run = 'EGP'
fileIndex = '502' #time in seconds

if len(fileIndex) == 3:
    filenumber = '0000'+str(fileIndex)
if len(fileIndex) == 4:
    filenumber = '000'+str(fileIndex)

f = pt.vlsvfile.VlsvReader(get_vlsvfile_fullpath(run, int(fileIndex)))
#f =pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.{}.vlsv".format(filenumber)) # vlsvReaderObject , from reading the .vlsv file
#root_dir = '/wrk-vakka/users/emilirin/plots/' # where you want to save the plot
#root_dir = '/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/new_plots/'
root_dir = '/wrk-vakka/users/horakons/carrington/plots/{}_validation_paper/'.format(run)

threshold = 1      # beta_star
delta = 0.05       # consider data satisfying beta_star = threshold +/- delta
                   # no use for delta anymore but I kept it here in case!!


#paper plot
#mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = True)  # magnetopause position at time t
#2023 proposal plot

#mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = True, plot_shue = False, scatter = False, plot_GEO = True, streamline_color = 'lime', streamline_density = 1.5, streamlinethick = 2, plotvar = 'vg_pe_egp')  # magnetopause position at time t
mp_nose_x = fit_magnetopause(f, run = run, root_dir = root_dir, fileIndex = fileIndex, threshold = threshold, delta = delta, plot = True, plot_shue = False, scatter = False, plot_GEO = True, streamline_color = 'lime', streamline_density = 1.5, streamlinethick = 2, plotvar = 'vg_j_jwc_egp')  # magnetopause position at time t


