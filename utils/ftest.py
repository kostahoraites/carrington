import pytools as pt

def f(vlsvfile = '/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.0001001.vlsv'):
   return pt.vlsvfile.VlsvReader(vlsvfile)
