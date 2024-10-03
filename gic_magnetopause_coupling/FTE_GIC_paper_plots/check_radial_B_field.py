#sanity check: ionosphere solver assumes dipole B-field is radial, but how radial is it?

import numpy as np
theta = 10. * np.pi / 180
r_hat = np.array([np.sin(theta), np.cos(theta)])
m = np.array([0, 1])
r1 = 3 * r_hat * (np.sum(m * r_hat) ) - m
r1_hat = r1 / np.linalg.norm(r1)
print('dipole magnetic field angle wrt radial vector, in degrees')
print( np.arccos(np.dot(r1_hat, r_hat)) * 180. / np.pi ) # Answer: ~5 degrees

