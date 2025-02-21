import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from myutils import timeseries, get_vlsvfile_fullpath
import ftest as ft

plt.rcParams.update({'font.size': 12})

# 1. estimate frequencies of different magnetopause processes, and relate these to the FTE-GIC time scaleu

# 1a) Freeman et al's magnetopause oscillations:

run = 'FHA'
dct = pd.read_csv('/wrk-vakka/users/horakons/carrington/magnetopause/ALL_CODES/mp_nose_x_new_FHA.csv')
start =1001
stop = 1496   #1496, 1612
t0 = 501 # initial time in the simulation

n_sw = 1e6 # SI, FHA parameter
v_F = 7.5e5  # SI, FHA parameter
m_p = 1.67e-27
R_E = 6.371e6
s = 1    # from Desai (2021) eq. 1, needed to match with Freeman

R_F = np.average(dct['x'][start-t0:]) * R_E
C_I = np.average(dct['C_I'][start-t0:])   # magnetosheath mass density integrated over sheath (units kg / m^2)
c = C_I * s / (n_sw * m_p * R_F) # dimensionless inertia term

# Horaites et al., page 11
b = v_F / (c * R_F)
omega = b * np.sqrt(6 * c - 1)
period = (2 * np.pi) / omega

print('Magnetopause oscillation (a la Freeman) period: {} seconds'.format(period))


# FFT of magnetopause oscillations
from scipy.fft import fft, fftfreq

import numpy as np


x_standoff_RE = dct['x'][start-t0:]
t = dct['t'][start-t0:]
ps = np.abs(np.fft.fft(x_standoff_RE))**2  # power spectrum

N = x_standoff_RE.size
time_step = 1.   # seconds
freqs = np.fft.fftfreq(N, time_step)

# sample spacing

import matplotlib.pyplot as plt

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(t, x_standoff_RE)
ax1.set_title(r'Magnetopause oscillations, $R(t)$')
ax1.set_xlabel(r'Time $[s]$')
ax1.set_ylabel(r'Subsolar standoff distance [$R_E$]')
ax2.plot(freqs, ps)
ax2.set_title('Power Spectral Density')
ax2.set_xlabel(r'Frequency f $[Hz]$')
ax2.set_ylabel(r'PSD $[R_E**2$/Hz]$')
ax2.set_yscale('log')

plt.tight_layout()
plt.savefig('/wrk-vakka/users/horakons/carrington/plots/{}/mp_io_coupling/mpause_standoff_fft.png'.format(run))
plt.close()


# 1b) Mirror mode waves:

#LOAD DATA
#start =1001
#stop = 1612   #1496, 1612

f_dummy = ft.f(get_vlsvfile_fullpath(run, start))

x_sc = 11.63; y_sc = 0; z_sc = 0.86
xvec_sc = np.array([x_sc, y_sc, z_sc]) * R_E

dct_ts = timeseries('FHA', ['proton/vg_rho', 'vg_b_vol.magnitude'], start, stop, x = xvec_sc)

'''
vsc_cellid = f_dummy.get_cellid(np.array([x_sc, y_sc, z_sc]) * 6.371e6)  # magnetopause at ~10.3 R_E
cellids = f_dummy.read_variable('cellid')
ind = np.where((cellids == vsc_cellid))[0][0]   # array index of the virtual spacecraft
'''
t = np.arange(start, stop+1)

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.set_size_inches(10, 4)

#vg_rho = dct_ts['proton/vg_rho'][:, ind] * 1e-6  # cm^-3
vg_rho = dct_ts['proton/vg_rho'] * 1e-6  # cm^-3
vg_rho_ave = np.average(vg_rho)
ax1.plot(t, vg_rho, color = 'blue', label = r'n [cm$^{-3}$]')
#ax1.plot(t, (vg_rho - vg_rho_ave) / vg_rho_ave, color = 'blue', label = r'$\Delta n/n$')
#bmag = dct_ts['vg_b_vol.magnitude'][:, ind]* 1e9      # nT
ax1.set_xlabel('time [s]')
ax1.set_ylabel(r'n [cm$^{-3}$]', color = 'blue')
delta_y = np.max([np.nanmax(vg_rho)-vg_rho_ave, vg_rho_ave - np.nanmin(vg_rho)])  # >0
ax1.set_ylim([vg_rho_ave - delta_y, vg_rho_ave + delta_y])
ax1b=ax1.twinx()
bmag = dct_ts['vg_b_vol.magnitude'] * 1e9      # nT
bmag_ave = np.average(bmag)
delta_y = np.max([np.nanmax(bmag)-bmag_ave, bmag_ave - np.nanmin(bmag)])  # >0
ax1b.set_ylim([bmag_ave - delta_y, bmag_ave + delta_y])
ax1b.plot(t, bmag, color = 'red', label = r'B [nT]')
ax1b.set_ylabel(r'B [nT]', color = 'red')
#ax1.plot(t, (bmag - bmag_ave) / bmag_ave, color = 'red', label = r'$\Delta B/B$')
#ax1.set_ylim([0, 40])
#ax1.set_ylim([-1, 1])
#ax1.legend()
ax1.set_title(r'Mirror modes at [{}, {}, {}]'.format(x_sc, y_sc, z_sc) + r' $R_E$')


ps1 = np.abs(np.fft.fft(vg_rho-vg_rho_ave))**2  # (mean-subtracted) power spectrum
ps2 = np.abs(np.fft.fft(bmag-vg_rho_ave))**2  # (mean-subtracted) power spectrum

N = vg_rho.size
time_step = 1.   # seconds
freqs = np.fft.fftfreq(N, time_step)
indfreqs = np.where(freqs > 0)[0] 

imax1 = np.argmax(ps1[indfreqs])
imax2 = np.argmax(ps2[indfreqs])

freqmax1 = freqs[indfreqs[imax1]] 
freqmax2 = freqs[indfreqs[imax2]] 

print('Mirror modes--- ')
print('Density spectrum peak: {} Hz'.format(freqmax1))
print('Magnetic spectrum peak: {} Hz'.format(freqmax2))

ax2.plot(freqs[indfreqs], ps1[indfreqs], color = 'blue', label = 'PSD(n-<n>)')
ax2.set_title('Power Spectral Density')
ax2.set_xlabel('Frequency f [Hz]')
ax2.set_ylabel(r'PSD(n-<n>) [cm$^{-6}$ Hz$^{-1}$]', color = 'blue')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.axvline(x = freqs[indfreqs[imax1]], color = 'blue', label = '{:.3f} Hz'.format(freqs[indfreqs[imax1]]), linestyle = '--')
ax2b = ax2.twinx() 
ax2b.plot(freqs[indfreqs], ps2[indfreqs], color = 'red', label = r'PSD(B-<B>) [nT$^2$ Hz$^{-1}$]')
ax2b.set_ylabel(r'PSD(B-<B>) [nT$^{2}$ Hz$^{-1}$]', color = 'red')
ax2b.axvline(x = freqs[indfreqs[imax2]], color = 'red', label = '{:.3f} Hz'.format(freqs[indfreqs[imax2]]), linestyle = '--')
ax2b.set_yscale('log')
#ax2.legend()
plt.tight_layout()

plt.savefig('/wrk-vakka/users/horakons/carrington/plots/{}/mp_io_coupling/mirror_mode_x{}_y{}_z{}.png'.format(run, x_sc, y_sc, z_sc))
plt.close()



