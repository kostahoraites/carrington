# example call: python coords_grid.py -nx 70 -nz 100 -xlim -15 25 -zlim -25 25

import numpy as np
import argparse

#Parse command line arguments:

parser = argparse.ArgumentParser()
parser.add_argument('-xlim', nargs='*', default=[0,0], help="xmin, xmax" )
parser.add_argument('-ylim', nargs='*', default=[0,0], help="ymin, ymax" )
parser.add_argument('-zlim', nargs='*', default=[0,0], help="zmin, zmax" )
parser.add_argument('-nx', default=1)
parser.add_argument('-ny', default=1)
parser.add_argument('-nz', default=1)

args = parser.parse_args()

xlim = np.array(args.xlim).astype(float)
ylim = np.array(args.ylim).astype(float)
zlim = np.array(args.zlim).astype(float)


# make a size [nx, ny, nz] grid
nx = int(args.nx)
ny = int(args.ny)
nz = int(args.nz)


# set corner point (coordinates in AU)

x0 = np.array([xlim[0], ylim[0], zlim[0]])

dx = np.array( [ (xlim[1]-xlim[0])/nx, 0, 0] )
dy = np.array( [ 0, (ylim[1]-ylim[0])/ny, 0] )
dz = np.array( [ 0, 0, (zlim[1]-zlim[0])/nz] )


output_file = open('txt_files/coords_grid.txt', 'w')

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            x= x0 + (i*dx) + (j*dy) + (k*dz)
            output_file.write(str(x).replace('[','').replace(']','') + '\n')


#for i in range(nr1):
#    for j in range(nr2):
#        x= x0 + (i*dr1) + (j*dr2)
#        output_file.write(str(x).replace('[','').replace(']','') + '\n')


output_file.close()


