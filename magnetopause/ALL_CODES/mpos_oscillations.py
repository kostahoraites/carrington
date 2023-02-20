import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt


#matplotlib.use("TkAgg")
fig, ax = plt.subplots(figsize=(6,4), dpi=240)

f = 2.44
c_D = 1.2
Q = 8.22e22 #Am^2
mu_0 = 4*np.pi*1e-7 #A/m^2
n = 4e6 #1/m^3
M = 1.67262192e-27 #kg (proton mass)
u = 7.5e5 #m/s
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_E = 6.371e6 #meters
R_01 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters
#R_01 = 10.25*R_E #meters


#t_shock = 857 # hard-code this for this crib sheet
t_shock = 860 # hard-code this for this crib sheet
#t_shock = 787 # hard-code this for this crib sheet

#STANDOFF DISTANCE (r) DATA:

t_r = []

for line in open('time_r_mp.txt', 'r'): #when step = 0.5
#for line in open('time_r_mp2.txt', 'r'): #when step = 1
    l = [i for i in line.split()]
    r = float(l[0])
    t = float(l[1])
    t_r.append([t, r])

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

f = 1.7
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_03 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters

#R_03 = 10.25 * R_E 

start = t_shock
stop = t_list[-1]
dt = 0.001
t3 = np.arange(start, stop, dt)
nt = t3.size + 1

x3 = np.zeros(nt)
v = np.zeros(nt)

#initial conditions:

x3[0] = R_03
v[0] = -0.013*R_E
#v[0] = -0.018*R_E

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
ax.plot(t3, x3/R_E, 'tab:cyan', label=r'Nonlinear, $f = 1.7$')

#NONLINEAR ^6 VS ^2:
ax.plot(t3, x3/R_E, 'tab:cyan', label=r'Nonlinear, $R^{-6}$ ($f = 1.7$)')
#ax.plot(t4, x4/R_E, 'tab:orange', label=r'Nonlinear, $R^{-2}$ ($f = 1.7$)')

#LINEAR:

#ax.plot(t_data, linear2/R_E, 'tab:red', linestyle='-.', alpha=0.5, label=r'Linear, $f = 2$')
#ax.plot(t_data, linear3/R_E, 'tab:cyan', linestyle='-.', alpha=0.5, label=r'Linear, $f = 1.7$')


#TITLES, AXES:
plt.grid(linewidth=0.5, alpha=0.5)
plt.title('Magnetopause Standoff Distance')
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


#MODEL ^6, f = 1.7:

f = 1.7
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_tdep = R_03   # (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters


start = t_shock
stop = t_list[-1]
dt = 0.001
t_tdep = np.arange(start, stop, dt)
nt = t_tdep.size + 1

x_tdep = np.zeros(nt)
v = np.zeros(nt)

#initial conditions:

x_tdep[0] = R_tdep
#v[0] = -0.013*R_E   # estimated from first 20 sec
#c_D = 1.2
v[0] = -0.018*R_E    # estimated from first 40 sec
#v[0] = 1e-3*R_E
#c_D = 2.2

# solve for x(t), v(t) using Euler's method

t0=80.  # sec
#n_nF = n_pulse(t_tdep-t_tdep[0], t0, n, type = 'tanh', tstart = 0) / n

t0_n = 911.37
t0_c = 960.98

delta = 35  # added this term to explain get pressure timing right

n_nF = (1.46e-9 * np.tanh((t_tdep-894. )/44.) +2.35e-9) / (1.46e-9 + 2.35e-9 )
#n_nF = (1.46e-9 * np.tanh((t_tdep-t0_n )/66.45) +2.35e-9) / (1.46e-9 + 2.35e-9 )


#c_D_tdep = c_D * n_pulse(t_tdep-t_tdep[0], t0*1.5, n, type = 'tanh') / n

toffset = 145 #  0, 75, 145
#c_D_tdep = 0.8 * np.tanh((t_tdep-t_tdep[0]-toffset)/95.)+1.32

c_D_tdep =  (0.48 * np.tanh((t_tdep-t0_c)/116.14)+0.87)

dt_pd_c = 85   # crossing time [sec] from 20 RE to magnetopause 
for i, t_i in enumerate(t_tdep):
    # dx/dt = v(t)
    x_tdep[i+1] = x_tdep[i] + (dt * v[i])
    # dv/dt =  -s/(c*R_D) * ((u+v)^2 - u^2*R_D^2/x^2)
    #v[i+1] = v[i] - dt/(c_D*R_D) * ( n_nF[i] * (u+v[i])**2 - u**2*R_D**6/x_tdep[i]**6)
    #v[i+1] = v[i] - dt/(c_D_tdep[i]*R_D) * ( n_nF[i] * (u+v[i])**2 - u**2*R_D**6/x_tdep[i]**6)
    #v[i+1] = v[i] - dt/(c_D * R_D) * ( (u+v[i])**2 - (1. / n_nF[i]) *u**2*R_D**6/x_tdep[i]**6)
    #if i > dt_pd_c:
    #    factor = n_nF[i-dt_pd_c]
    #else:
    #    factor = 0.25    
    #v[i+1] = v[i] - dt/(c_D * R_D * factor) * ( n_nF[i] * (u+v[i])**2 - u**2*R_D**6/x_tdep[i]**6)
    v[i+1] = v[i] - dt/(c_D_tdep[i] * R_D) * ( n_nF[i] * (u+v[i])**2 - u**2*R_D**6/x_tdep[i]**6)



t_tdep = np.array(list(t_tdep) + [t_tdep[-1]+dt])



#STANDOFF DISTANCE DATA:
ax.plot(t_list[:shock_index], r_list[:shock_index]/R_E, color='k', linestyle=':')
ax.plot(t_data, r_data/R_E, 'k', label='Data')

#R_0:

ax.plot([621, t_shock], [R_tdep/R_E, R_tdep/R_E], color='tab:orange', linestyle=':')

#NONLINEAR:

ax.plot(t_tdep, x_tdep/R_E, 'tab:orange', label=r'Nonlinear, $f = 1.7$, time-dep.')

#NONLINEAR ^6 VS ^2:
ax.plot(t_tdep, x_tdep/R_E, 'tab:orange', label=r'Nonlinear, $R^{-6}$ ($f = 1.7$), time-dep.')

#LINEAR:

#ax.plot(t_data, linear3/R_E, 'tab:orange', linestyle='-.', alpha=0.5, label=r'Linear, $f = 1.7$')

#-------------------------------------------------------------------------------

#plt.show()

plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/magnetopause_oscillations_fig_b.png')
plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/magnetopause_oscillations_fig_b.pdf')

