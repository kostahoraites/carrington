import pytools as pt
import numpy as np
from myutils import mkdir_path


save_dir = '/wrk-vakka/users/horakons/carrington/data/particle_tracer/vlsv/'
f_dummy = pt.vlsvfile.VlsvReader( '/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/bulk/dense_cold_hall1e5_afterRestart374/bulk1.0001000.vlsv' )

# ExB drift ~ 10^(E_exp - B_exp) m/s
Ex = 0; Ey = -1; Ez = 0
E_exp = -3

Bx = 0; By = 0; Bz = 1
B_exp = -8


E0 = np.array([Ex, Ey, Ez]) * 10.**E_exp  # V/m
B0 = np.array([Bx, By, Bz]) * 10.**B_exp  # T

fg_e = f_dummy.read_variable('fg_e')*0 + E0
fg_b = f_dummy.read_variable('fg_b')*0 + B0


vg_b_vol = f_dummy.read_variable('vg_b_vol')*0 + B0
vg_e_vol = vg_b_vol*0 + E0

B_mag = np.repeat(np.linalg.norm(vg_b_vol, axis = -1)[:, None], 3, axis = 1)

# set proton velocity to be exactly ExB drift
vg_v = np.cross(vg_e_vol, vg_b_vol) / B_mag**2
vg_rho = f_dummy.read_variable('proton/vg_rho')*0 + 1e6  # as EGI SW density

cellid = f_dummy.read_variable("CellID")

# write to file
filename_vlsv = save_dir + 'E_{}_{}_{}_e{}_B_{}_{}_{}_e{}.vlsv'.format(Ex, Ey, Ez, E_exp, Bx, By, Bz, B_exp)
mkdir_path(filename_vlsv)
writer = pt.vlsvfile.VlsvWriter(f_dummy, filename_vlsv)
#writer.write(pos,'ig_r','SpatialGrid','ionosphere')
writer.write(fg_e.reshape([fg_e.size//3, 3]),'fg_e','VARIABLE','fsgrid')  # fg_ grid arrays need to be reshaped to dimensions [ncells,3] 
writer.write(fg_b.reshape([fg_b.size//3, 3]),'fg_b','VARIABLE','fsgrid')
writer.write(vg_b_vol,'vg_b_vol','VARIABLE','SpatialGrid')
writer.write(vg_v,'proton/vg_v','VARIABLE','SpatialGrid')
writer.write(vg_rho,'proton/vg_rho','VARIABLE','SpatialGrid')
writer.write(cellid,'CellID','VARIABLE','SpatialGrid')
