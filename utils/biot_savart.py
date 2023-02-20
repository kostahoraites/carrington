import pytools as pt
from pyCalculations.intpol_points import vlsv_intpol_points
import pyCalculations.fieldtracer as fieldtracer
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from static_field_tracer_3d_alt import static_field_tracer_3d_alt
from myutils import *    #e.g. this imports get_vlsvfile_fullpath, mkdir_path, cartesian_to_spherical, spherical_to_cartesian, numcurl3d, numjacobian3d
import os, sys
import argparse
import warnings
from copy import deepcopy
import scipy
import pandas
from time import time


global R_EARTH
R_EARTH = 6.371e6            #check what value is used in simulations
#global CELLSIZE
global ROOT_DIR
ROOT_DIR = '/wrk-vakka/users/horakons/carrington/plots/'
#ROOT_DIR = '/wrk-vakka/users/horakons/carrington/test_plots/'
global mu_0
mu_0 = 4e-7 * np.pi




