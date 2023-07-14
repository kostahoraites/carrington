#python 3 required

import pandas as pd
from pyamps import AMPS

 # EGI and EGL parameters (note n not used)

for Bz, n, run in zip([-5., -10.], [1., 4.], ['Control', 'Pulse']):
    m = AMPS(750, # Solar wind velocity in km/s
             0, # IMF By (GSM) in nT
             Bz, # IMF Bz (GSM) in nT,
             0, # dipole tilt angle in degrees
             178.4) # F107_index in September 2001, see Tulasi Ram et al. 2009
    Ju = m.get_upward_current()   # calculate field-aligned currents on a pre-defined grid
    mlat, mlt = m.scalargrid    # Ju will be evaluated at the following coords:
    df = pd.DataFrame({'Jd': -Ju.flatten(), 'mlat': mlat.flatten(),  'mlt': mlt.flatten()})  # Jd = -Ju, Downward current
    df.to_csv('pyamp_{}.csv'.format(run), index=False, sep = ' ')
