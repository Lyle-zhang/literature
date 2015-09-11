"""
Compare 1-D transient heat conduction model to Papadikis2010a for d=350um
"""

import numpy as np
import matplotlib.pyplot as py
from transhc import hc
from kinetics import kn2

# Parameters
#------------------------------------------------------------------------------

rhow = 700      # density of wood, kg/m^3
d = 0.035e-2    # biomass particle diameter, m
cpw = 1500      # biomass specific heat capacity, J/kg*K
cpc = 1100      # char specific heat capacity, J/kg*K
kw = 0.105      # biomass thermal conductivity, W/m*K
kc = 0.071      # char thermal conductivity, W/m*K
h = 900         # heat transfer coefficient, W/m^2*K
Ti = 300        # initial particle temp, K
Tinf = 773      # ambient temp, K
H = 255000     # heat of reaction, J/kg

# Shape factor, time, and node (radius point) vectors
#------------------------------------------------------------------------------

b = 2       # run model as a cylinder (b = 1) or as a sphere (b = 2)
nt = 2000                       # number of time steps
tmax = 1                        # max time, s
dt = tmax/nt                    # time step, s
t = np.arange(0, tmax+dt, dt)   # time vector

nr = 19     # number or radius steps
r = d/2     # radius of particle, m
dr = r/nr   # radius step, delta r
m = nr+1    # nodes from center m=0 to surface m=steps+1

# Temperature and Density arrays, Mass Fraction vector
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
pt = np.zeros((len(t), m))      # create array for tar density

pw[0] = rhow                 # initial wood density at all nodes

# mass fraction vector
# columns = average mass fraction of entire solid at a time step
Ys = np.zeros(len(t))   # create row vector for mass fraction
Ys[:] = 1               # initial mass fraction, Ys=1 for all wood

# Initial thermal properties 
#------------------------------------------------------------------------------

Yw = pw[0]/rhow     # wood fraction, Yw=1 all wood, Yw=0 all char

cpbar = Yw*cpw + (1-Yw)*cpc     # effective heat capacity
kbar = Yw*kw + (1-Yw)*kc        # effective thermal conductivity
pbar = pw[0] + pc[0]            # effective density

g = np.ones(m)*(1e-10)  # assume initial heat generation is negligible

# Solve system of equations [A]{T}={C} where T = A\C for each time step
#------------------------------------------------------------------------------

for i in range(1, nt+1):
    
    # heat conduction
    T[i] = hc(m, dr, b, dt, h, Tinf, g, T, i, r, pbar, cpbar, kbar)
    
    # kinetic reactions
    pw[i], pc[i], pg[i], pt[i], g = kn2(T, pw, pc, pg, pt, dt, i, H)
    
    # update mass fraction vector
    Yw = pw[i] / rhow
    cpbar = Yw*cpw + (1-Yw)*cpc
    kbar = Yw*kw + (1-Yw)*kc
    pbar = pw[i] + pc[i]
    Ys[i] = np.mean(pbar) / rhow
    
Tavg = [np.mean(row) for row in T]  # average temperature for entire particle

# Plot
#------------------------------------------------------------------------------

# grab data from csv file
tc350, Tc350 = np.loadtxt('Fig7_cent350.csv', delimiter=',', unpack=True)
ts350, Ts350 = np.loadtxt('Fig7_surf350.csv', delimiter=',', unpack=True)

# plot model vs data
py.close('all')

py.figure(1)
py.plot(t, T[:, 0], '-b', label='center')
py.plot(t, T[:, m-1], '-r', label='surface')
py.plot(t, Tavg, '-m', label='avg')
py.plot(tc350-1, Tc350, 'ob', label='center')
py.plot(ts350-1, Ts350, 'or', label='surface')
py.axhline(Tinf, c='k', ls='--', label='ambient')
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')
py.title(r'Temperatures for d={:.0f}$\mu m$, h={}$W/m^2K$'.format(d*10**6, h))
py.legend(loc='best', numpoints=1)
py.ylim(ymin=Ti-20)

py.grid()
py.show()