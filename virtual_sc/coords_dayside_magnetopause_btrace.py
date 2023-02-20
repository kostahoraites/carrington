
import numpy as np
from static_field_tracer3d import static_field_tracer3d
import ftest
from myutils import get_vlsvfile_fullpath
from carrington import inner_radius

R_EARTH = 6.371e6

run='egl'
fileIndex=800
#x0=[-6.5*R_EARTH, 0.8*R_EARTH, 1.9*R_EARTH]
x0=[9*R_EARTH, 0*R_EARTH, 0*R_EARTH]   # dayside magnetopause, roughly

#run='egi'
#fileIndex=1489
#x0=[-7*R_EARTH, 0*R_EARTH, 3.2*R_EARTH]

file = get_vlsvfile_fullpath(run, fileIndex)
#file = '/wrk-vakka/group/spacephysics/vlasiator/3D/EGL/bulk/bulk1.egl.0001760.vlsv'
f = ftest.f(file)


max_iterations = 1500
dx = 0.01 * R_EARTH
ncoordsave = 150    # total number of coordinates in the final file
dstep_write = int(2*max_iterations / ncoordsave)   # number of iterations between writes

#x_b = static_field_tracer3d( f, x0, max_iterations, dx, direction='+-', bvar='fg_b' )
x_b = static_field_tracer3d( f, x0, max_iterations, dx, direction='+-', bvar='fg_b' )


###
# cut out anything inside inner radius:

to_pop = []
min_radius = inner_radius(f)
for i, x in enumerate(x_b):
    r = np.sqrt(x[0]**2 + x[1]**2 + x[2]**2)
    if r/R_EARTH < min_radius:
        to_pop = to_pop + [i]
        print(r/R_EARTH)

to_pop.reverse()  # pop the largest indices first

for i in to_pop:
    x_b.pop(i)
###

output_file = open('txt_files/coords_dayside_magnetopause_btrace_{}.txt'.format(run), 'w')
i = 0
for x in x_b:
    if i == dstep_write:
        output_file.write(str(x/R_EARTH).replace('[','').replace(']','') + '\n')
        i=0
    else:
        i+=1
output_file.close()


