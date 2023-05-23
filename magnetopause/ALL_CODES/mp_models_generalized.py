import numpy as np
import matplotlib
import matplotlib.pyplot as plt


#matplotlib.use("TkAgg")
fig, ax = plt.subplots(figsize=(6,4), dpi=240)


f = 2.44
c_D = 1.2
Q = 8e22 #Am^2   # Vlasiator dipole
#Q = 8.22e22 #Am^2      # Earth's real dipole
mu_0 = 4*np.pi*1e-7 #A/m^2
n = 4e6 #1/m^3
M = 1.67262192e-27 #kg (proton mass)
u = 7.5e5 #m/s
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_E = 6.371e6 #meters
R_01 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters
#R_01 = 10.25*R_E #meters



#v_0 = 0.    # estimated from first 40 sec
#v_0 = -0.013*R_E   # estimated from first 20 sec
#v_0 = -0.018*R_E    # estimated from first 40 sec
v_0 = -0.018 * R_E

#t_shock = 857 # hard-code this for this crib sheet
t_shock = 857 # hard-code this for this crib sheet
#t_shock = 787 # hard-code this for this crib sheet

#STANDOFF DISTANCE (r) DATA:

t_r = []

for line in open('time_r_mp.txt', 'r'): #when step = 0.5
#for line in open('time_r_mp2.txt', 'r'): #when step = 1
    l = [i for i in line.split()]
    r = float(l[0])
    t_temp = float(l[1])
    t_r.append([t_temp, r])

t_r.sort() #sorted list of [t, r] elements

t_list = []
r_list = []
R_E = 6.371e6 #meters

for i in t_r:
    t_list.append(i[0]) #list of t values (621–1760 s)
    r_list.append(i[1]*R_E) #list of r values in meters

t_list = np.array(t_list)
r_list = np.array(r_list)

#shock_index = np.where(r_list == np.max(r_list[200:]))[0][-1] #236
shock_index = int(t_shock - 621)
r_shock = r_list[shock_index]
t_shock = t_list[shock_index]
r_data = r_list[shock_index:] #list of r values after shock
t_data = t_list[shock_index:] #list of t values after shock (857–1760 s)

t_data_0 = []

for i in t_data:
    i_0 = i-t_data[0]
    t_data_0.append(i_0)

t_data_0 = np.array(t_data_0) #list of t values starting from zero (0–903 s)
r_data = np.array(r_data)


#-------------------------------------------------------------------------------
#MODEL ^6, f = 1.7:

#f = 1.7
f = 1.75
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_03 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters

#R_03 = 10.25 * R_E 

start = t_shock
stop = t_list[-1]
dt = 1e-3   # 1e-6
t3 = np.arange(start, stop, dt)
nt = t3.size + 1

x3 = np.zeros(nt)
v = np.zeros(nt)

#initial conditions:

x3[0] = R_03
v[0] = v_0

# solve for x(t), v(t) using Euler's method

for i, t_i in enumerate(t3):
    # dx/dt = v(t)
    x3[i+1] = x3[i] + (dt * v[i])
    # dv/dt =  -s/(c*R_D) * ((u+v)^2 - u^2*R_D^2/x^2)
    v[i+1] = v[i] - dt/(c_D*R_D) * ((u+v[i])**2 - u**2*R_D**6/x3[i]**6)

t3 = np.array(list(t3) + [t3[-1]+dt])


#----------------------------------------------------------------------------


#STANDOFF DISTANCE DATA:
ax.plot(t_list[:shock_index], r_list[:shock_index]/R_E, color='k', linestyle=':')
ax.plot(t_data, r_data/R_E, 'k', label='Data')

#R_0:

ax.plot([621, t_shock], [R_03/R_E, R_03/R_E], color='tab:cyan', linestyle=':')
#ax.plot([621, t_shock], [R_04/R_E, R_04/R_E], color='tab:orange', linestyle=':')

#NONLINEAR:
#ax.plot(t, x/R_E, 'goldenrod', label=r'Nonlinear, $f = 2.44$')
ax.plot(t3, x3/R_E, 'tab:cyan', label=r'Eq. 3')

#NONLINEAR ^6 VS ^2:
#ax.plot(t3, x3/R_E, 'tab:cyan', label=r'Nonlinear, $R^{-6}$ ($f = 1.7$)')
#ax.plot(t4, x4/R_E, 'tab:orange', label=r'Nonlinear, $R^{-2}$ ($f = 1.7$)')

#LINEAR:

#ax.plot(t_data, linear2/R_E, 'tab:red', linestyle='-.', alpha=0.5, label=r'Linear, $f = 2$')
#ax.plot(t_data, linear3/R_E, 'tab:cyan', linestyle='-.', alpha=0.5, label=r'Linear, $f = 1.7$')


#TITLES, AXES:
plt.grid(linewidth=0.5, alpha=0.5)
#plt.title('Magnetopause Standoff Distance')
plt.xlabel('Time [s]')
ax.set_xticks([800, 857, 1000, 1200, 1400, 1600])
ax.set_xticklabels(['800', r'$t_{0}$', '1000', '1200', '1400', '1600'])
plt.ylabel(r'Standoff Distance R [$R_{E}$]')
plt.xlim(t_list[0], t_list[-1])



#-------------------------------------------------------------------------------
#NEW STUFF

import numpy as np

def n_pulse(t, t0, nF, type ='tanh', tstart=0):
    # t=time, t0=ramp timescale, nF = final density [SI]
    if type=='tanh':
        return nF * (0.75 * np.tanh((t-tstart) / t0) + 0.25)
    if type=='sigmoid':
        return nF * (0.75 * np.exp((t-tstart) / t0) / (1. + np.exp(t / t0)) + 0.25)


#MODEL ^6, f = 51.7:

#f = 1.7
f = 1.75
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_tdep = R_03   # (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters


start = t_shock
stop = t_list[-1]
t = np.arange(start, stop, dt)
nt = t.size + 1

x_tdep = np.zeros(nt)
v = np.zeros(nt)

#initial conditions:
import copy

x_tdep[0] = R_tdep
x_tdep_fit = copy.deepcopy(x_tdep)
v[0] = v_0

# solve for x(t), v(t) using Euler's method

###
# APPROACH 1: USE eta(t), c(t) FROM FITS
###

##eta = (1.46e-9 * np.tanh((t-t0_n )/(66.45)) +2.35e-9) / (1.46e-9 + 2.35e-9 )
eta_fit = 0.38 * np.tanh((t-943.52 )/87.21) +0.62

c_fit =  0.48 * np.tanh((t-961.01)/116.14) +0.87

for i, t_i in enumerate(t):
    # dx/dt = v(t)
    x_tdep_fit[i+1] = x_tdep_fit[i] + (dt * v[i])
    v[i+1] = v[i] - dt/(c_fit[i] * R_D) * ( eta_fit[i] * (u+v[i])**2 - u**2*R_D**6/x_tdep_fit[i]**6)



###
# APPROACH 2: INTERPOLATE eta(t), c(t)
###

import pandas as pd
from scipy import interpolate

df = pd.read_csv('/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/mp_nose_x_new_EGL.csv')

n_F = 4e6  # [SI]
m_p = 1.67e-27
RE = 6371000.
RF = 8.22 * RE
c_data = df['C_I'] / (m_p * n_F * RF)
eta_data = df['n'] / 3.96e6
#eta_data = df['n'] / n_F


t_cross = 35   # crossing time [sec] from 20 RE to magnetopause 

#smooth
#from scipy.signal import savgol_filter
#c_data = savgol_filter(c_data, 51, 2)
#eta_data = savgol_filter(eta_data, 51, 2)

f_c = interpolate.interp1d(df['t'], c_data, fill_value='extrapolate', bounds_error = False)
f_eta = interpolate.interp1d(df['t']+t_cross, eta_data, fill_value='extrapolate', bounds_error = False)
#f_u = interpolate.interp1d(df['t']+t_cross, df['v'], fill_value='extrapolate', bounds_error = False)

c = f_c(t)   # use interpolation function returned by `interp1d`
eta = f_eta(t)   # use interpolation function returned by `interp1d`
#u_i = f_u(t)

###

for i, t_i in enumerate(t):
    # dx/dt = v(t)
    x_tdep[i+1] = x_tdep[i] + (dt * v[i])
    v[i+1] = v[i] - dt/(c[i] * R_D) * ( eta[i] * (u+v[i])**2 - u**2*R_D**6/x_tdep[i]**6)
    #v[i+1] = v[i] - dt/(c[i] * R_D) * ( eta[i] * (u_i[i]+v[i])**2 - u**2*R_D**6/x_tdep[i]**6)



t = np.array(list(t) + [t[-1]+dt])



#STANDOFF DISTANCE DATA:
ax.plot(t_list[:shock_index], r_list[:shock_index]/R_E, color='k', linestyle=':')
#ax.plot(t_data, r_data/R_E, 'k', label='Data')

#R_0:

ax.plot([621, t_shock], [R_tdep/R_E, R_tdep/R_E], color='tab:orange', linestyle=':')

#NONLINEAR (original model):

#ax.plot(t, x_tdep/R_E, 'tab:orange', label=r'Eq. 3')

#NONLINEAR (generalized model):
ax.plot(t, x_tdep_fit/R_E, 'tab:green', label=r'Eq. 6: $\eta$, $c$ sigmoid')
ax.plot(t, x_tdep/R_E, 'tab:orange', label=r'Eq. 6: $\eta$, $c$ interpolated')

#LINEAR:

#ax.plot(t_data, linear3/R_E, 'tab:orange', linestyle='-.', alpha=0.5, label=r'Linear, $f = 1.7$')

#-------------------------------------------------------------------------------

ax.legend()

#ax.set_yticks([7.5, 8, 8.5, 9, 9.5, 10,10.5,11,11.5,12])

plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/original_vs_generalized.pdf')

