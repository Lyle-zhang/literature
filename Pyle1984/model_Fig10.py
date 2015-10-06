"""
Compare 1-D transient heat conduction model to Figure 10 in Pyle 1984
"""

import numpy as np
import matplotlib.pyplot as py
from transhc import hc
from kinetics import kn
    
# Parameters
#------------------------------------------------------------------------------

rhow = 450      # density of wood, kg/m^3
d = 0.006       # biomass particle diameter, m
h = 55          # heat transfer coefficient, W/m^2*K
Ti = 303        # initial particle temp, K
Tinf = 780      # ambient temp, K
H = -100000     # heat of reaction, J/kg

# Shape factor, time, and node (radius point) vectors
#------------------------------------------------------------------------------

b = 1           # run model as a cylinder (b = 1) or as a sphere (b = 2)
nt = 2000                       # number of time steps
tmax = 150                      # max time, s
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector

nr = 19     # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

# Temperature and Density arrays, Conversion vector
#------------------------------------------------------------------------------

# temperture array
# rows = time step, columns = node points from center to surface node
T = np.zeros((len(t), m))   # create array to store temperatures
T[0] = Ti                   # initial temperature at all nodes

# density array
# rows = time step, columns = node points from center to surface
pw = np.zeros((len(t), m))      # create array for wood density
pc = np.zeros((len(t), m))      # create array for char density
pg = np.zeros((len(t), m))      # create array for gas density

pw[0] = rhow                    # initial wood density at all nodes

# mass fraction array
B = np.ones((len(t), m))
C1 = np.zeros((len(t), m))
C2 = np.zeros((len(t), m))

# conversion vector
# columns = average conversion of entire solid at a time step
Ys = np.zeros(len(t))   # create row vector for conversion, Ys=0 for all wood

# Initial thermal properties 
#------------------------------------------------------------------------------

Yw = pw[0]/rhow     # wood fraction, Yw=1 all wood, Yw=0 all char

cpw = 1112.0 + 4.85 * (T[0] - 273.15)   # wood heat capacity, J/(kg*K) 
kw = 0.13 + (3e-4) * (T[0] - 273.15)    # wood thermal conductivity, W/(m*K)
cpc = 1003.2 + 2.09 * (T[0] - 273.15)   # char heat capacity, J/(kg*K)
kc = 0.08 - (1e-4) * (T[0] - 273.15)    # char thermal conductivity, W/(m*K)

cpbar = Yw*cpw + (1-Yw)*cpc             # effective heat capacity
kbar = Yw*kw + (1-Yw)*kc                # effective thermal conductivity
pbar = pw[0] + pc[0]                    # effective density

g = np.ones(m)*(1e-10)  # assume initial heat generation is negligible

# Solve system of equations [A]{T}={C} where T = A\C for each time step
#------------------------------------------------------------------------------

for i in range(1, nt+1):
    
    # heat conduction
    T[i] = hc(m, dr, b, dt, h, Tinf, g, T, i, r, pbar, cpbar, kbar)
    
    # kinetic reactions
    B[i], C1[i], C2[i], g = kn(T, B, C1, C2, rhow, dt, i, H)
    
    # update thermal properties
    cpw = 1112.0 + 4.85 * (T[i] - 273.15)
    kw = 0.13 + (3e-4) * (T[i] - 273.15)
    cpc = 1003.2 + 2.09 * (T[i] - 273.15)
    kc = 0.08 - (1e-4) * (T[i] - 273.15)
    
    # update wood and char density
    pw[i] = B[i]*rhow
    pc[i] = (C1[i]+C2[i])*rhow
    
    # update mass fraction vector
    Yw = pw[i] / (pw[i] + pc[i])
    cpbar = Yw*cpw + (1-Yw)*cpc
    kbar = Yw*kw + (1-Yw)*kc
    pbar = pw[i] + pc[i]
    Ys[i] = 1 - np.mean(pw[i])/rhow
    
Tavg = [np.mean(row) for row in T]  # average temperature for entire particle

# Plot
#------------------------------------------------------------------------------

# preferences for plots
py.rcParams['lines.linewidth'] = 2
py.rcParams['axes.grid'] = True

# grab experiment data from csv file
t1, Temp = np.loadtxt('Fig10Tcenter.csv', delimiter=',', unpack=True)
t2, Mass = np.loadtxt('Fig10conv.csv', delimiter=',', unpack=True)

# plot model vs data
py.close('all')

py.figure(1)
py.plot(t/60, T[:, 0], '-b', label='center')
py.plot(t1, Temp, 'og', mec='g', label='expmt')
py.axhline(Tinf, c='k', ls='--', label='ambient')
py.ylim(ymin=Ti-20)
py.legend(loc='best', numpoints=1)
py.xlabel('Time (min)')
py.ylabel('Temperature (K)')
#py.title(r'Temperatures for d={:.0f} $mm$, h={} $W/m^2K$'.format(d*1000, h))

# remove top, right spines and tick marks
ax = py.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

py.figure(2)
py.plot(t/60, Ys, '-b', label='model')
py.plot(t2, Mass, 'og', mec='g', label='expmt')
py.axhline(1, c='k', ls='--')
py.ylim([0, 1.1])
py.legend(loc='best', numpoints=1)
py.xlabel('Time (min)')
py.ylabel('Conversion (-)')

ax = py.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

py.show()
