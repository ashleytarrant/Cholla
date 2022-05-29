#!/usr/bin/env python3
# Example python plotting script for the 1D Sod Shock Tube test

import h5py
import numpy as np
import matplotlib
matplotlib.rcParams['mathtext.default']='regular'
matplotlib.rcParams['xtick.direction']='in'
matplotlib.rcParams['ytick.direction']='in'
matplotlib.rcParams['xtick.top']=True
matplotlib.rcParams['ytick.right']=True
import matplotlib.pyplot as plt

mp = 1.67e-24
kb = 1.38e-16

dnamein='./'
dnameout='./'

DE = 1 # dual energy flag - 1 if the test was run with dual energy
i = 1 # output file number

f = h5py.File(dnamein+str(i)+'.h5', 'r')
head = f.attrs
nx = head['dims'][0]
dx = head['dx'][0]
gamma = head['gamma'][0]
dunit = head['density_unit']
vunit = head['velocity_unit']
punit = head['energy_unit']
d  = np.array(f['density']) # mass density
mx = np.array(f['momentum_x']) # x-momentum
my = np.array(f['momentum_y']) # y-momentum
mz = np.array(f['momentum_z']) # z-momentum
E  = np.array(f['Energy']) # total energy density
vx = mx/d
vy = my/d
vz = mz/d
if DE:
  e  = np.array(f['GasEnergy'])
  p  = e*(gamma-1.0)
  ge = e/d
else: 
  p  = (E - 0.5*d*(vx*vx + vy*vy + vz*vz)) * (gamma - 1.0)
  ge  = p/d/(gamma - 1.0)
n = d*dunit / (0.6*mp)
vx = vx*vunit / 1e5
p = p*punit
T = p / (n*kb)
p = p/kb
log_p = np.log10(p)
log_T = np.log10(T)

print(np.mean(n), np.mean(vx), np.mean(p))
xplot = np.arange(nx)*dx

fig = plt.figure(figsize=(6,6))
ax1 = plt.axes([0.12, 0.6, 0.35, 0.35])
plt.axis([0, nx*dx, 0.001, 0.01])
ax1.plot(xplot, n[:,64,64], 'o', markersize=2, color='black')
plt.xlabel('x [kpc]')
plt.ylabel(r'$n$ [cm$^{-3}$]')
ax2 = plt.axes([0.6, 0.6, 0.35, 0.35])
plt.axis([0, nx*dx, 0.0, 1100])
ax2.plot(xplot, vx[:,64,64], 'o', markersize=2, color='black')
plt.xlabel('x [kpc]')
plt.ylabel(r'$v$ [km / s]')
ax3 = plt.axes([0.12, 0.1, 0.35, 0.35])
plt.axis([0, nx*dx, 4, 5])
ax3.plot(xplot, log_p[:,64,64], 'o', markersize=2, color='black')
plt.xlabel('x [kpc]')
plt.ylabel(r'log(P/k) [K cm$^{-3}$]')
ax4 = plt.axes([0.6, 0.1, 0.35, 0.35])
plt.axis([0, nx*dx, 6.7, 7.0])
ax4.plot(xplot, log_T[:,64,64], 'o', markersize=2, color='black')
plt.xlabel('x [kpc]')
plt.ylabel(r'log(T) [K]')

plt.savefig(dnameout+str(i)+".png", dpi=300);
plt.close(fig)
