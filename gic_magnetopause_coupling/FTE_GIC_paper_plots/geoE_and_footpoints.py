import pytools as pt


fileIndex = 1300
filename = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1_sidecars/geoelectric_field/ionosphere_gic_FHA_{}.vlsv'.format(str(fileIndex).zfill(7))
pt.plot.plot_ionosphere(filename=filename, outputdir="./", var="ig_E_north", lin=True, wmark="NE", viewdir=1)

'''
filename = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.0001100.vlsv'
pt.plot.plot_ionosphere(filename=filename, outputdir="./", var="ig_fac", lin=True, wmark="NE", viewdir=1)
'''
