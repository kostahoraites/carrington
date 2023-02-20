import numpy as np
import warnings
from scipy import interpolate


def static_field_tracer_3d( vlsvReader, coord_list, max_iterations, dx, direction='+', fg = None ):
    ''' static_field_tracer_3d() integrates along the (static) field-grid vector field to calculate a final position. 
        Code uses forward Euler method to conduct the tracing.
        Based on Analysator's static_field_tracer()

        :Inputs:
         param vlsvReader:      A vlsvReader object (~an open .vlsv file)
         param coord_list:      a list of 3-element array-like initial coordinates [ [x1,y1,z1], [x2,y2,z2], ... ]
                                if considering just a single starting point, the code accepts a 3-element array-like object [x1,y1,z1]
         param max_iterations:  The maximum amount of iterations before the algorithm stops
         param dx:              One iteration step length (ex. dx=R_EARTH/50)
         keyword direction:     '+' or '-' or '+-' Follow field in the plus direction or minus direction
         keyword fg:            The field grid variable to be traced (either a string or numpy array)
                                options include:
                                    fg = some string
                                        ex. fg='fg_b': B-field, fg='fg_e': E-field
                                        static_field_tracer_3d() will load the appropriate variable via the vlsvReader object
                                    fg = None
                                        same as inputting fg='fg_b' above
                                    fg = some field-grid ("fg") array.          dimensions [dimx,dimy,dimz,3]
                                        ex. fg = vlsvobj.read_variable('fg_b')
                                        field grid data is already loaded externally using read_variable() method (see vlsvreader.py).
                                   
        :returns:                  points_traced --- Traced coordinates (a list of lists of 3-element coordinate arrays)
                                   ex. points_traced[2][5][1]: at 3rd tracing step [2], the 6th point [5], y-coordinate [1]
                                      note: Can convert output to a 3D numpy array if desired, with np.array(points_traced)

        EXAMPLE:            vlsvobj = pytools.vlsvfile.VlsvReader(vlsvfile) 
                            fg_b = vlsvobj.read_variable('fg_b')
                            static_field_tracer_3d( vlsvobj, [[5e7,0,0], [0,0,5e7]], 10, 1e5, direction='+', fg = fg_b )

    '''

    # standardize input (a list of 3-element arrays/lists)
    if type(coord_list[0]) not in [list, np.ndarray]:
        coord_list = [coord_list]

    # Read face-centered field grid variable (denoted 'fg_*' in .vlsv files):
    if type(fg) == str:
        fg = vlsvReader.read_variable(fg)
    elif (fg is None):
        fg = vlsvReader.read_variable('fg_b')        # e.g., for EGL run: fg.shape = (1024, 736, 736, 3)
    else:
        #   fg is already an ndarray
        warnings.warn("User defined tracing variable detected. fg, volumetric variable results may not work as intended, use face-values instead.")
    
    # Recursion (trace in both directions and sum the result)
    if direction == '+-':
        backward = static_field_tracer_3d(vlsvReader, coord_list, max_iterations, dx, direction='-', fg=fg)
        backward.reverse()
        forward = static_field_tracer_3d(vlsvReader, coord_list, max_iterations, dx, direction='+', fg=fg)
        return backward + forward[1:]

    # Create x, y, and z coordinates:
    xsize = fg.shape[0]
    ysize = fg.shape[1]
    zsize = fg.shape[2]
    xmin = vlsvReader.read_parameter('xmin')
    xmax = vlsvReader.read_parameter('xmax')
    ymin = vlsvReader.read_parameter('ymin')
    ymax = vlsvReader.read_parameter('ymax')
    zmin = vlsvReader.read_parameter('zmin')
    zmax = vlsvReader.read_parameter('zmax')
    sizes = np.array([xsize, ysize, zsize])
    maxs = np.array([xmax, ymax, zmax])
    mins = np.array([xmin, ymin, zmin])
    dcell = (maxs - mins)/(sizes.astype('float'))
    x = np.arange(mins[0], maxs[0], dcell[0]) + 0.5*dcell[0]
    y = np.arange(mins[1], maxs[1], dcell[1]) + 0.5*dcell[1]
    z = np.arange(mins[2], maxs[2], dcell[2]) + 0.5*dcell[2]
    coordinates = np.array([x,y,z], dtype=object)

    # Create grid interpolation of vector field (V)
    interpolator_face_V_0 = interpolate.RegularGridInterpolator((x-0.5*dcell[0], y, z), fg[:,:,:,0], bounds_error = False, fill_value = np.nan)
    interpolator_face_V_1 = interpolate.RegularGridInterpolator((x, y-0.5*dcell[1], z), fg[:,:,:,1], bounds_error = False, fill_value = np.nan)
    interpolator_face_V_2 = interpolate.RegularGridInterpolator((x, y, z-0.5*dcell[2]), fg[:,:,:,2], bounds_error = False, fill_value = np.nan)
    interpolators = [interpolator_face_V_0, interpolator_face_V_1, interpolator_face_V_2]

    # Trace vector field lines
    if direction == '-':
        multiplier = -1
    else:
        multiplier = 1
    points_traced = [np.array(coord_list)]              # iteratively append traced trajectories to this list
    points = points_traced[0]
    N = len(list(coord_list))
    V_unit = np.zeros([3, N])
    for i in range(max_iterations):
        V_unit[0, :] = interpolators[0](points)
        V_unit[1, :] = interpolators[1](points)
        V_unit[2, :] = interpolators[2](points)
        V_mag = np.linalg.norm(V_unit, axis=(0))
        V_unit = V_unit / V_mag[np.newaxis,:]
        new_points = points + multiplier*V_unit.T * dx
        points = new_points
        points_traced.append( list(points) )             # list of lists of 3-element arrays
    return points_traced


