import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
#-------------------------------------------------------------------------------
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

shock_index = np.where(r_list == np.max(r_list[200:]))[0][-1] #236
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
#MODEL ^6, f = 2.44:

def func(x, t):
    return 0                           # no forcing
    #return np.exp(-np.abs(x)/(1+t))

# d^2x/ dt^2 + 1/(c_D*R_D) * ((u+dx/dt)^2 - u^2*R_D^2/x^2) = func(x,t)

# rewrite the above equation as 2 equations:
# v(t) = dx/dt
# dv/dt + 1/(c_D*R_D) * ((u+v)^2 - u^2*R_D^2/x^2) = func(x, t)

# Pick some coefficients

f = 2.44
c_D = 1.2
Q = 8.22e22 #Am^2
mu_0 = 4*np.pi*1e-7 #A/m^2
n = 4e6 #1/m^3
M = 1.67262192e-27 #kg (proton mass)
u = 7.5e5 #m/s
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_01 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters
R_E = 6.371e6 #meters

# define time, x, v arrays

start = t_shock
stop = t_list[-1]
dt = 1e-3
t = np.arange(start, stop, dt)
nt = t.size + 1

x = np.zeros(nt)
v = np.zeros(nt)


# use Scipy's ODE integrator (LSODA method, detects stiffness)
from scipy.integrate import odeint

# let y = [R, dR/dt]
def desai(y, t, c, v_F, R_F):
    R, dRdt = y
    dydt = [dRdt, - 1./(c*R_F) * ((v_F+dRdt)**2 - v_F**2*R_F**6/R**6)]
    return dydt

#initial conditions:

x[0] = R_01
v[0] = -0.018*R_E    # linear trend for first 40 seconds
#v[0] = -0.013*R_E   # linear trend for first 20 seconds

# solve for x(t), v(t) using Euler's method

for i, t_i in enumerate(t):
    # dx/dt = v(t)
    x[i+1] = x[i] + (dt * v[i])
    # dv/dt =  -s/(c*R_D) * ((u+v)^2 - u^2*R_D^2/x^2)
    v[i+1] = v[i] - dt/(c_D*R_D) * ((u+v[i])**2 - u**2*R_D**6/x[i]**6)


t = np.array(list(t) + [t[-1]+dt])


# use Scipy's ODE integrator
#sol = odeint(desai, [R_01, -0.018*R_E], t, args=(c_D, u, R_D), rtol = 1e-12, atol = 1e-12)
#x = sol[:, 0]     # v = sol[:,1] 

#-------------------------------------------------------------------------------
#MODEL ^6, f = 2:

f = 2
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_02 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters

start = t_shock
stop = t_list[-1]
t2 = np.arange(start, stop, dt)
nt = t2.size + 1

x2 = np.zeros(nt)
v = np.zeros(nt)



#initial conditions:

x2[0] = R_02
v[0] = -0.018*R_E    # linear trend for first 40 seconds
#v[0] = -0.013*R_E   # linear trend for first 20 seconds

# solve for x(t), v(t) using Euler's method

for i, t_i in enumerate(t2):
    # dx/dt = v(t)
    x2[i+1] = x2[i] + (dt * v[i])
    # dv/dt =  -s/(c*R_D) * ((u+v)^2 - u^2*R_D^2/x^2)
    v[i+1] = v[i] - dt/(c_D*R_D) * ((u+v[i])**2 - u**2*R_D**6/x2[i]**6)

t2 = np.array(list(t2) + [t2[-1]+dt])




# use Scipy's ODE integrator
#sol = odeint(desai, [R_02, -0.018*R_E], t2, args=(c_D, u, R_D), rtol = 1e-12, atol = 1e-12)
#x2 = sol[:, 0]     # v = sol[:,1] 


#-------------------------------------------------------------------------------
#MODEL ^6, f = 1.7:

f = 1.7
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_03 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters

start = t_shock
stop = t_list[-1]
t3 = np.arange(start, stop, dt)
nt = t3.size + 1

x3 = np.zeros(nt)
v = np.zeros(nt)

#initial conditions:

x3[0] = R_03
v[0] = -0.018*R_E    # linear trend for first 40 seconds
#v[0] = -0.013*R_E   # linear trend for first 20 seconds

# solve for x(t), v(t) using Euler's method

for i, t_i in enumerate(t3):
    # dx/dt = v(t)
    x3[i+1] = x3[i] + (dt * v[i])
    # dv/dt =  -s/(c*R_D) * ((u+v)^2 - u^2*R_D^2/x^2)
    v[i+1] = v[i] - dt/(c_D*R_D) * ((u+v[i])**2 - u**2*R_D**6/x3[i]**6)

t3 = np.array(list(t3) + [t3[-1]+dt])

# use Scipy's ODE integrator
#sol = odeint(desai, [R_03, -0.018*R_E], t3, args=(c_D, u, R_D), rtol = 1e-12, atol = 1e-12)
#x3 = sol[:, 0]     # v = sol[:,1] 

#-------------------------------------------------------------------------------
#MODEL ^2, f = 1.7:

f = 1.7
R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.) #meters
R_04 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.) #meters

start = t_shock
stop = t_list[-1]
t4 = np.arange(start, stop, dt)
nt = t4.size + 1

x4 = np.zeros(nt)
v = np.zeros(nt)

#initial conditions:

x4[0] = R_04
v[0] = -0.018*R_E    # linear trend for first 40 seconds
#v[0] = -0.013*R_E   # linear trend for first 20 seconds
#v[0] = -0.02*R_E

# solve for x(t), v(t) using Euler's method

#c_D = 1.2
#eta = -1

for i, t_i in enumerate(t4):
    # dx/dt = v(t)
    x4[i+1] = x4[i] + (dt * v[i])
    # dv/dt =  -s/(c*R_D) * ((u+v)^2 - u^2*R_D^2/x^2)
    v[i+1] = v[i] - dt/(c_D*R_D) * ((u+v[i])**2 - u**2*R_D**2/x4[i]**2)
    #v[i+1] = v[i] - dt/(c_D*R_D) * ((u+v[i])**2 - u**2*R_D**(-2*eta)/x4[i]**(-2*eta))

t4 = np.array(list(t4) + [t4[-1]+dt])


# use Scipy's ODE integrator
#sol = odeint(desai, [R_04, -0.018*R_E], t4, args=(c_D, u, R_D), rtol = 1e-12, atol = 1e-12)
#x4 = sol[:, 0]     # v = sol[:,1] 

#------------------------------------------------------------------------------
#LINEAR MODEL:

f_values = [2.44, 2, 1.7]
for f in f_values:
    R_D = (f**2*mu_0*Q**2/(32*np.pi**2*n*M*u**2))**(1/6.)
    R_0 = (f**2*mu_0*Q**2/(32*np.pi**2*1e6*M*u**2))**(1/6.)
    ampl = R_0 - R_D
    tau = R_D/u
    K = 1.2
    b = 1/(K*tau) # ~ 0.012
    phi = 0
    if f == 2.44:
        linear1 = R_D + ampl*np.exp(-b*t_data_0)*np.cos(b*(6*K-1)**0.5*t_data_0 + phi)
    if f == 2:
        linear2 = R_D + ampl*np.exp(-b*t_data_0)*np.cos(b*(6*K-1)**0.5*t_data_0 + phi)
    if f == 1.7:
        linear3 = R_D + ampl*np.exp(-b*t_data_0)*np.cos(b*(6*K-1)**0.5*t_data_0 + phi)

#-----------------------------------------------------------------------------
#SHUE 98:

D_0 = (1e6*M*u**2)*1e9
D_inf = (n*M*u**2)*1e9 #nPa
B_0 = -5 #nT
B_inf = -10 #nT
alpha_0 = (0.58 - 0.007 * B_0) * (1 + 0.024 * np.log(D_0))
r_0 = (10.22 + 1.29 * np.tanh(0.184 * (B_0 + 8.14))) * D_0 ** (-1/6.6)
alpha_inf = (0.58 - 0.007 * B_inf) * (1 + 0.024 * np.log(D_inf))
r_inf = (10.22 + 1.29 * np.tanh(0.184 * (B_inf + 8.14))) * D_inf ** (-1/6.6)

#----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6,4), dpi=240)

#STANDOFF DISTANCE DATA:
ax.plot(t_list[:shock_index], r_list[:shock_index]/R_E, color='k', linestyle=':')
ax.plot(t_data, r_data/R_E, 'k', label='Data')

#R_0:
ax.plot([621, t_shock], [R_01/R_E, R_01/R_E], color='goldenrod', linestyle=':')
ax.plot([621, t_shock], [R_02/R_E, R_02/R_E], color='tab:red', linestyle=':')
ax.plot([621, t_shock], [R_03/R_E, R_03/R_E], color='tab:cyan', linestyle=':')
#ax.plot([621, t_shock], [R_04/R_E, R_04/R_E], color='tab:orange', linestyle=':')

#NONLINEAR:
ax.plot(t, x/R_E, 'goldenrod', label=r'Nonlinear, $f = 2.44$')
ax.plot(t2, x2/R_E, 'tab:red', label=r'Nonlinear, $f = 2$')
ax.plot(t3, x3/R_E, 'tab:cyan', label=r'Nonlinear, $f = 1.7$')

#NONLINEAR ^6 VS ^2:
#ax.plot(t3, x3/R_E, 'tab:cyan', label=r'Nonlinear, $R^{-6}$ ($f = 1.7$)')
#ax.plot(t4, x4/R_E, 'tab:orange', label=r'Nonlinear, $R^{-2}$ ($f = 1.7$)')

#LINEAR:
ax.plot(t_data, linear1/R_E, 'goldenrod', linestyle='-.', alpha=0.5, label=r'Linear, $f = 2.44$')
ax.plot(t_data, linear2/R_E, 'tab:red', linestyle='-.', alpha=0.5, label=r'Linear, $f = 2$')
ax.plot(t_data, linear3/R_E, 'tab:cyan', linestyle='-.', alpha=0.5, label=r'Linear, $f = 1.7$')

#SHUE:
ax.plot([621, t_shock], [r_0, r_0], color='limegreen', linestyle='--', linewidth=1.5, label='Shue (1998)')
ax.plot([1400, 1760], [r_inf, r_inf], color='limegreen', linestyle='--', linewidth=1.5)

#TITLES, AXES:
plt.grid(linewidth=0.5, alpha=0.5)
plt.title('Magnetopause Standoff Distance')
plt.xlabel('Time [s]')
ax.set_xticks([800, 857, 1000, 1200, 1400, 1600])
ax.set_yticks([7.5, 8, 8.5, 9, 9.5, 10,10.5,11,11.5,12])
ax.set_xticklabels(['800', r'$t_{0}$', '1000', '1200', '1400', '1600'])
plt.ylabel(r'Standoff Distance R [$R_{E}$]')
plt.xlim(t_list[0], t_list[-1])
#ax.legend()

plt.legend( loc = "upper right")
plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/nonlinear_vs_linear.pdf')
#plt.savefig('all_models.pdf')

plt.close()


