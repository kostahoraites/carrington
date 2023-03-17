from myutils import *
from read_cdf import read_cdf
import numpy as np
import matplotlib.pyplot as plt
from mpos_oscillations import mpos_oscillations


data = read_cdf('omni_hro_1min_20150301_v01.cdf')

Epoch = data['Epoch']
BX_GSE = data['BX_GSE']
BY_GSE = data['BY_GSE']
BZ_GSE = data['BZ_GSE']
flow_speed = data['flow_speed']
Vx = data['Vx']
Vy = data['Vy']
Vz = data['Vz']
proton_density = data['proton_density']
T = data['T']
Pressure = data['Pressure']
Mach_num = data['Mach_num']
Mgs_mach_num = data['Mgs_mach_num']
AU_INDEX = data['AU_INDEX']
SYM_H = data['SYM_H']

#t_start = t_double('2015-03-17/04:30:00', epoch=True)
#t_stop = t_double('2015-03-17/05:00:00', epoch=True)
t_start = t_double('2015-03-17/00:00:00', epoch=True)
t_stop = t_double('2015-03-18/00:00:00', epoch=True)

mp = 1.67e-27
R_E = 6371000
P_d = mp * (proton_density * 1e6) * (flow_speed * 1e3)**2. * 1e9       # nPa

#ind, = np.where((Epoch >= t_start) & (Epoch <= t_stop))
ind, = np.where((Epoch >= t_start) & (Epoch <= t_stop) & (BZ_GSE < 9e3) & (flow_speed < 9e4) & (proton_density < 999) & (P_d < 1e6) & (SYM_H < 9e3))
t = Epoch[ind]
t_minutes = (t -t_start) / (1000 * 60)


names = [r'$B_Z$ [GSE]', 'v [km/sec]', r'$n [cm^{-3}]$', r'$P_d$ [nPa]', 'SYM-H [nT]']
variables = [BZ_GSE[ind], flow_speed[ind], proton_density[ind], P_d[ind], SYM_H[ind]] 
dummy = [9999, 9e4, 999, 1e6, 9999]

nvar = len(variables)
fig, axes = plt.subplots(nvar, 1, figsize= (12,7))



for i in range(nvar-1):
    ind_plot, = np.where(variables[i] < dummy[i])
    axes[i].plot(t_minutes[ind_plot], variables[i][ind_plot])
    axes[i].set_ylabel(names[i])
    axes[i].set_xlabel('')
    for xlabel_i in axes[i].get_xticklabels():
        xlabel_i.set_visible(False)


# SYM-H
axes[-1].plot(t_minutes[ind_plot], variables[-1][ind_plot], color = "C1")    #last row
axes[-1].set_xlabel('minutes after ' + t_string(t_start, epoch=True)[0:19])
axes[-1].set_ylabel(names[-1])


#magnetopause
t_in = (t - t[0]) / 1000
n = proton_density[ind] * 1e6
v = flow_speed[ind] * 1e3
c = 1
R, t_R = mpos_oscillations(t_in, n, v, 1, dt = 1e-3, t_cross = 0, R_0 = None, dRdt_0 = 0, f=2.44, dip_mom = 8.22e22)

ax2=axes[-1].twinx()
# make a plot with different y-axis using second axis object
ax2.plot(t_R/60, R/R_E, color = "C2", label = r'Theoretical $R(t)$')  #,marker="o")
ax2.set_ylabel(r'$R(t) [R_E]$',color="C2",fontsize=18)
ax2.set_ylim([7,10])





plt.savefig('st_patricks_storm.png')
plt.close()








