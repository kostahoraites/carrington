# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import numpy as np
import scipy as sp
import pytools as pt
import warnings
from scipy import interpolate

def dynamic_field_tracer( vlsvReader_list, x0, max_iterations, dx):
   ''' Field tracer in a dynamic time frame

      :param vlsvReader_list:              List of vlsv readers
      :param x0:                           The starting point for the streamlines
      
   '''
   dt = vlsvReader_list[1].read_parameter('t') - vlsvReader_list[0].read_parameter('t')
   # Loop through vlsvreaders:
   iterations = 0
   for vlsvReader in vlsvReader_list:
      v = vlsvReader.read_interpolated_variable('v', x0)
      stream_plus = static_field_tracer( vlsvReader, x0, max_iterations, dx, direction='+' )
      stream_minus = static_field_tracer( vlsvReader, x0, max_iterations, dx, direction='-' )
      stream = stream_minus[::-1] + stream_plus # Minus reversed
      pt.miscellaneous.write_vtk_file("test" + str(iterations) + ".vtk", stream)
      x0 = x0 + v*dt
      iterations = iterations + 1

def static_field_tracer( vlsvReader, x0, max_iterations, dx, direction='+', bvar='B' ):
   ''' Field tracer in a static frame

       :param vlsvReader:         An open vlsv file
       :param x:                  Starting point for the field trace
       :param max_iterations:     The maximum amount of iteractions before the algorithm stops
       :param dx:                 One iteration step length
       :param direction:          '+' or '-' or '+-' Follow field in the plus direction or minus direction
       :param bvar:               String, variable name to trace [default 'B']
       :returns:                  List of coordinates
   '''

   if(bvar is not 'B'):
     warnings.warn("User defined tracing variable detected. fg, volumetric variable results may not work as intended, use face-values instead.")

   if direction == '+-':
     backward = static_field_tracer(vlsvReader, x0, max_iterations, dx, direction='-', bvar=bvar)
     backward.reverse()
     forward = static_field_tracer(vlsvReader, x0, max_iterations, dx, direction='+', bvar=bvar)
     return backward + forward

   f = vlsvReader
   # Read cellids in order to sort variables
   cellids = vlsvReader.read_variable("CellID")
   xsize = f.read_parameter("xcells_ini")
   ysize = f.read_parameter("ycells_ini")
   zsize = f.read_parameter("zcells_ini")
   xmin = f.read_parameter('xmin')
   xmax = f.read_parameter('xmax')
   ymin = f.read_parameter('ymin')
   ymax = f.read_parameter('ymax')
   zmin = f.read_parameter('zmin')
   zmax = f.read_parameter('zmax')

   sizes = np.array([xsize, ysize, zsize])
   maxs = np.array([xmax, ymax, zmax])
   mins = np.array([xmin, ymin, zmin])
   dcell = (maxs - mins)/(sizes.astype('float'))

   # Pick only two coordinate directions to operate in
   if xsize <= 1:
      indices = [2,1]
   elif ysize <= 1:
      indices = [2,0]
   elif zsize <= 1:
      indices = [1,0]
   else:
      if bvar is 'B':
         bvar = 'fg_b'
      return static_field_tracer3d( vlsvReader, x0, max_iterations, dx, direction=direction, bvar=bvar )
      
   # Read face_B:
   face_B = f.read_variable(bvar)
   face_Bx = face_B[:,0]
   face_By = face_B[:,1]
   face_Bz = face_B[:,2]

   face_Bx = face_Bx[cellids.argsort()].reshape(sizes[indices])
   face_By = face_By[cellids.argsort()].reshape(sizes[indices])
   face_Bz = face_Bz[cellids.argsort()].reshape(sizes[indices])

   face_B = np.array([face_Bx, face_By, face_Bz])

   # Create x, y, and z coordinates:
   x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
   y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
   z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
   coordinates = np.array([x,y,z])
   # Debug:
   if( len(x) != sizes[0] ):
      print("SIZE WRONG: " + str(len(x)) + " " + str(sizes[0]))

   # Create grid interpolation
   interpolator_face_B_0 = interpolate.RectBivariateSpline(coordinates[indices[0]] - 0.5*dcell[indices[0]], coordinates[indices[1]], face_B[indices[0]], kx=2, ky=2, s=0)
   interpolator_face_B_1 = interpolate.RectBivariateSpline(coordinates[indices[0]], coordinates[indices[1]] - 0.5*dcell[indices[1]], face_B[indices[1]], kx=2, ky=2, s=0)
   interpolators = [interpolator_face_B_0, interpolator_face_B_1]#, interpolator_face_B_2]

   #######################################################
   if direction == '-':
      multiplier = -1
   else:
      multiplier = 1

   points = [x0]
   for i in range(max_iterations):
      previous_point = points[-1]
      B_unit = np.zeros(3)
      B_unit[indices[0]] = interpolators[0](previous_point[indices[0]], previous_point[indices[1]])
      B_unit[indices[1]] = interpolators[1](previous_point[indices[0]], previous_point[indices[1]])
      B_unit = B_unit / float(np.linalg.norm(B_unit))
      points.append( previous_point + multiplier*B_unit * dx )
   #######################################################

   return points


def static_field_tracer3d( vlsvReader, x0, max_iterations, dx, direction='+', bvar='fg_b' ):
   ''' Field tracer in a static frame
       :param vlsvReader:         An open vlsv file
       :param x:                  Starting point for the field trace, 3-element numpy array
       :param max_iterations:     The maximum amount of iteractions before the algorithm stops
       :param dx:                 One iteration step length
       :param direction:          '+' or '-' or '+-' Follow field in the plus direction or minus direction
       :param bvar:               String, variable name to trace [default 'fg_b'].
                                  OR an ndarray, i.e. obtained vlsvReader.reader_parameter() externally and passed into the function
       :returns:                  List of coordinates
   '''

   if(bvar is not 'fg_b'):
     warnings.warn("User defined tracing variable detected. fg, volumetric variable results may not work as intended, use face-values instead.")

   if direction == '+-':
     backward = static_field_tracer3d(vlsvReader, x0, max_iterations, dx, direction='-', bvar=bvar)
     backward.reverse()
     forward = static_field_tracer3d(vlsvReader, x0, max_iterations, dx, direction='+', bvar=bvar)
     return backward + forward

   f = vlsvReader
   # Read cellids in order to sort variables
   cellids = vlsvReader.read_variable("CellID")
   xmin = f.read_parameter('xmin')
   xmax = f.read_parameter('xmax')
   ymin = f.read_parameter('ymin')
   ymax = f.read_parameter('ymax')
   zmin = f.read_parameter('zmin')
   zmax = f.read_parameter('zmax')

   # Read face_B:
   if type(bvar) == np.ndarray:
      face_B = bvar
   else:
      face_B = f.read_variable(bvar)

   xsize = face_B.shape[0]
   ysize = face_B.shape[1]
   zsize = face_B.shape[2]

   sizes = np.array([xsize, ysize, zsize])
   maxs = np.array([xmax, ymax, zmax])
   mins = np.array([xmin, ymin, zmin])
   dcell = (maxs - mins)/(sizes.astype('float'))

   # Create x, y, and z coordinates:
   x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
   y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
   z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
   # Debug:
   if( len(x) != sizes[0] ):
      print("SIZE WRONG: " + str(len(x)) + " " + str(sizes[0]))

   # Create grid interpolation

   interpolator_face_B_0 = interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), face_B[:,:,:,0])
   interpolator_face_B_1 = interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), face_B[:,:,:,1])
   interpolator_face_B_2 = interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), face_B[:,:,:,2])
   interpolators = [interpolator_face_B_0, interpolator_face_B_1, interpolator_face_B_2]

   #######################################################
   if direction == '-':
      multiplier = -1
   else:
      multiplier = 1

   points = [x0]
   for i in range(max_iterations):
      previous_point = points[-1]
      B_unit = np.zeros(3)
      B_unit[0] = interpolators[0]([previous_point[0], previous_point[1], previous_point[2]])
      B_unit[1] = interpolators[1]([previous_point[0], previous_point[1], previous_point[2]])
      B_unit[2] = interpolators[2]([previous_point[0], previous_point[1], previous_point[2]])
      B_unit = B_unit / float(np.linalg.norm(B_unit))
      points.append( previous_point + multiplier*B_unit * dx )
   #######################################################

   return points






