import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd

import scipy.optimize

from read_vsc_data import *    # read_vsc_data reads txt file output from vlsvintpol.py into a dictionary

#matplotlib.use("TkAgg")

fig,ax = plt.subplots()

run = 'EGL'

#dct_vsc, nvars_vsc = read_vsc_data('/wrk-vakka/users/horakons/carrington/virtual_sc/txt_files/virtual_sc_data_vlasiator_comp_{}.txt'.format(run.lower()))
#dct_vsc, nvars_vsc = read_vsc_data('/wrk-vakka/users/horakons/carrington/virtual_sc/txt_files/virtual_sc_data_vlasiator_comp_{}_draft1.txt'.format(run.lower()))

#df = pd.read_csv('/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/mp_nose_x_new.csv')
df = pd.read_csv('/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/mp_nose_x_new_EGL.csv')


df.keys()
x = df['x']
B=df['B']
t = df['t']
C_I = df['C_I']
pdyn = df['P_d'] * 1e9 # [nPa], at magnetopause + 5RE
#pdyn = dct_vsc['proton/vg_pdyn'][:,0]  * 1e9    # [nPa], at 15 RE  (20 RE in draft1 file)
n = df['n'] * 1e-6   # [cm^-3]
v = df['v']


c=np.arange(0,x.size)

n_F = 4e6  # [SI]
m_p = 1.67e-27
RE = 6371000.
RF = 8.22 * RE
c = C_I / (m_p * n_F * RF)

# every item is a 2d array [x,y] where x indexes time, y indexes position

t_vsc = t
#t_vsc = dct_vsc['t'][:,0]
X_RE = x
#X_RE = dct_vsc['X_RE'][0,:]


ntime = t_vsc.size
#npos_vsc = X_RE.size
npos_vsc = x.size


####### DYNAMIC PRESSURE

#t_offset = 100            # data from x=R+5RE
t_offset = 35            # data from x=15
t_pdyn = t_vsc + t_offset

pdyn_max = np.max(pdyn)



def func(t, A, t_0, t_T, y_0):
    return y_0 + A*np.tanh( (t-t_0)/t_T )

# CURVE FITTING
#fit to dynamic pressure

p0 = [1, 800, 100, 1]
p = scipy.optimize.curve_fit(func, t_pdyn, pdyn, p0=p0)[0]
P_d_model = func(t_pdyn, p[0], p[1], p[2], p[3])
A, t_0, t_T, y_0 = p
####

tmodel = np.arange(1000)
#ax.plot(t_pdyn, pdyn, label=r'Vlasiator $P_d(t-??)$, x=15 RE', color ="C2")     #, marker='o')
#ax.plot(t_pdyn, P_d_model, label ='{:.2f} * tanh((t-{:.2f})/{:.2f}) + {:.2f}'.format(A, t_0, t_T, y_0 ), color = "C2", linestyle = '--')  # , marker = "o")
ax.set_xlabel('t [seconds]', fontsize = 18)
#ax.set_ylabel(r'$P_d$ [nPa]',color="C2",fontsize=18)
ax.set_title(r'$P_d(t-35)$ and $c(t)$, Pulse Run')
ax.set_ylim([0,4])
#ax.legend(loc = 'center right', fontsize = 11)



####### DENSITY

#t_offset = 100            # data from x=R+5RE
t_offset = 35            # data from x=15
t_n = t + t_offset




def func(t, A, t_0, t_T, y_0):
    return y_0 + A*np.tanh( (t-t_0)/t_T )

# CURVE FITTING
#fit to dynamic density

p0 = [1, 800, 100, 1]
p = scipy.optimize.curve_fit(func, t_n, n, p0=p0)[0]
n_model = func(t_n, p[0], p[1], p[2], p[3])
A, t_0, t_T, y_0 = p
####

tmodel = np.arange(1000)
ax.plot(t_n, n, label=r'Vlasiator $n_2(t-35)$', color = "C0")     #, marker='o')
ax.plot(t_n, n_model, label ='{:.2f} * tanh((t-{:.2f})/{:.2f}) + {:.2f}'.format(A, t_0, t_T, y_0 ), color = "C0", linestyle = '--')  # , marker = "o")
ax.set_xlabel('t [seconds]', fontsize = 18)
ax.set_ylabel(r'$n_2 [cm^{-3}]$',color="C0",fontsize=18)
ax.set_title(r'$n_2(t-35)$ and $c(t)$, Pulse Run', fontsize = 20)
ax.set_ylim([0,4])
ax.legend(loc = 'center right', fontsize = 11)

####### C_I




# create figure and axis objects with subplots()


#SECOND PLOT
# twin object for two different y-axis on the sample plot
ax2=ax.twinx()
# make a plot with different y-axis using second axis object
ax2.plot(t, c,color = "C1", label = r'Vlasiator $c(t)$')  #,marker="o")
ax2.set_ylabel(r'$c$ [unitless]',color="C1",fontsize=18)
ax2.set_ylim([0,1.6])

#### CURVE FITTING
#fit to inertia c(t)

p0 = [1., 800, 100, 1.]
p = scipy.optimize.curve_fit(func, t, c, p0=p0)[0]
c_model = func(t, p[0], p[1], p[2], p[3])

A, t_0, t_T, y_0 = p
####

ax2.plot(t, c_model, label = '{:.2f} * tanh((t-{:.2f})/{:.2f}) + {:.2f}'.format(A, t_0, t_T, y_0 ), color = "C1", linestyle = '--')  # , marker = "o")


# save the plot as a file
#fig.savefig('two_diffeREnt_y_axis_for_single_python_plot_with_twinx.jpg',
#            format='jpeg',
#            dpi=100,
#            bbox_inches='tight')

plt.axvline(x = 857, color = 'black', label = r'$t_0$ = 857 s', linestyle = ':')
ax2.legend(loc='lower right', fontsize = 11)

#plt.show()


plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/Pd_c.png')
plt.savefig('/wrk-vakka/users/horakons/carrington/plots/EGL_validation_paper/Pd_c.pdf')

