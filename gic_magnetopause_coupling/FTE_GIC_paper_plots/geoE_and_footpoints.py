import matplotlib.pyplot as plt
import matplotlib.projections as projections
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pytools as pt

#fig, ax = plt.subplots()

figsize = (6,5)

minlatitude = 60

fig = plt.figure(figsize=figsize,dpi=150)
ax_cartesian = fig.add_axes([0.1,0.1,0.9,0.9], xlim=(-(90-minlatitude),(90-minlatitude)), ylim=(-(90-minlatitude),(90-minlatitude)), aspect='equal')
#ax_polar = fig.add_axes([0.1,0.1,0.9,0.9], polar=True, frameon=False, ylim=(0, 90-minlatitude))
ax_polar = inset_axes(parent_axes=ax_cartesian, width="100%", height="100%", axes_class = projections.get_projection_class('polar'), borderpad=0)
ax_polar.set_frame_on(False)
ax_polar.set_aspect('equal')


fileIndex = 1100
filename = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/geoelectric_field/ionosphere_gic_FHA_{}.vlsv'.format(str(fileIndex).zfill(7))
pt.plot.plot_ionosphere(filename=filename, outputdir="./", lin = True, scale = 2.0, axes = ax_polar, vscale = 1e3, vmin = -0.1, vmax = 0.1,
                        outputfile='ig_E_north_t_{}'.format(fileIndex), var="ig_E_north", viewdir=1, cbtitle = r'$E_\lambda$ [V/km]')

plt.savefig('./ig_E_north_t_{}.png'.format(fileIndex))
plt.close()

filename = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.0001100.vlsv'
pt.plot.plot_ionosphere(filename=filename, outputdir="./", outputfile='ig_fac_t'.format(fileIndex), var="ig_fac", lin=True, viewdir=1)
