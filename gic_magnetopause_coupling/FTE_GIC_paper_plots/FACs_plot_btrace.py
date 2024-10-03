import pytools as pt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from myutils import get_vlsvfile_fullpath

R_E = 6.371e6

run = 'FHA'
fileIndex = 1165
ax = plt.gca()
bulkname = get_vlsvfile_fullpath(run, fileIndex)

cmap = matplotlib.cm.plasma
cmap.set_under('grey')

outputdir = '/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/FTE_GIC_paper_plots/'
pt.plot.plot_colormap3dslice(filename=bulkname,var='proton/vg_pressure', boxre=[0, 14, -14, 14], normal = 'y', run=run,
                            colormap='plasma',step=fileIndex,outputdir=outputdir, vmin=0,vmax=1e-9,lin = True,
                            outputfile='FACs_plot_NO_btrace.pdf',
                            Earth=1, streamlines='vg_b_vol', streamlinedensity=2.5, streamlinethick = 2, streamlinecolor = 'black',
                            cutpointre=0, axes = ax, scale=1.5, useimshow=True)


x=np.loadtxt('/wrk-vakka/users/horakons/carrington/gic_magnetopause_coupling/txt_files/fg_innerboundary_btrace_C2_FHA.txt')

ax.plot(x[:,0]/R_E, x[:,2]/R_E, color = 'green', linewidth=5)
plt.savefig('FACs_plot_btrace.png')



