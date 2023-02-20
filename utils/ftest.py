from carrington import *
import pytools as pt


def f(vlsvfile = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.0001760.vlsv"):
   return pt.vlsvfile.VlsvReader(vlsvfile)
