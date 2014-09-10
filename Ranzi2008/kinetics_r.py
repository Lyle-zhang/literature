# Kinetics from Ranzi2008 paper
# plots species, C, as function of time, t

# use Python 3 print function
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

#---- global parameters

cell = 0.5      # cellulose
hemi = 0.3      # hemicellulose
lig = 0.2       # lignin

R = 1.987       # gas constant, kcal/kmol*K
T = 773        # temperature, K

dt = 0.002                       # time step, s
tmax = 1                        # time max, s
t = np.arange(0, tmax+dt, dt)   # time range, s
#t = np.linspace(0, 1)
nt = len(t)                     # number of time steps, -

#---- kinetic scheme from Table 1: Cellulose

# A = pre-exponential factor, 1/s and E = activation energy, kcal/kmol
A1 = 8e13;      E1 = 46000
A2 = 1e9;       E2 = 30000
A3 = 4;         E3 = 10000
A4 = 8e7;       E4 = 32000

# reaction rates, Arrhenius equation K = A exp(-E/RT)
K1 = A1 * np.exp(-E1 / (R * T))   
K2 = A2 * np.exp(-E2 / (R * T))
K3 = A3 * T * np.exp(-E3 / (R * T))
K4 = A4 * np.exp(-E4 / (R * T))

#---- kinetic scheme from Table 2: Hemicellulose

A5 = 1e10;      E5 = 31000
A6 = 3e9;       E6 = 27000
A7 = 3;         E7 = 11000
A8 = 1e10;      E8 = 33000

K5 = A5 * np.exp(-E5 / (R * T))
K6 = A6 * np.exp(-E6 / (R * T))
K7 = A7 * T * np.exp(-E7 / (R * T))
K8 = A8 * np.exp(-E8 / (R * T))

#---- kinetic scheme from Table 3: Lignin

A9 = 4e15;      E9 = 48500
A10 = 2e13;     E10 = 37500
A11 = 1e9;      E11 = 25500
A12 = 5e6;      E12 = 31500
A13 = 1e13;     E13 = 49500
A14 = 1e5;      E14 = 20500
A15 = 8e1;      E15 = 12000
A16 = 1.2e9;    E16 = 30000

K9 = A9 * np.exp(-E9 / (R * T))
K10 = A10 * np.exp(-E10 / (R * T))
K11 = A11 * np.exp(-E11 / (R * T))
K12 = A12 * np.exp(-E12 / (R * T))
K13 = A13 * np.exp(-E13 / (R * T))
K14 = A14 * np.exp(-E14 / (R * T))
K15 = A15 * T * np.exp(-E15 / (R * T))
K16 = A16 * np.exp(-E16 / (R * T))

#---- reactions for Table 1: Cellulose

r1 = np.zeros(nt)
r1[0] = -(K1 + K4)*cell

CELL = np.zeros(nt)
CELL[0] = cell

CELLA = np.zeros(nt)
CELLA[0] = 0

for i in range(1, nt):
    r1[i] = -(K1 + K4)*CELL[i-1]
    
    CELL[i] = CELL[i-1] - (K1 + K4)*CELL[i-1]*dt
    CELLA[i] = CELLA[i-1] + K1*CELL[i]*dt

#---- plot T vs K for cellulose, hemicellulose, lignin

# testing
tt = np.linspace(0, 1)
c = cell*np.exp(-(K1+K4)*tt)

py.figure(1)
py.plot(tt, c)
py.title('test')
py.show()

# main plot
py.figure(2)
py.plot(t, CELL)
py.plot(t, CELLA)
py.title('main')
py.show()



#py.figure(1)
#py.plot(t, r1, label='r1')
#py.legend(loc='best', numpoints=1)
#py.xlabel('Time (s)')
#py.ylabel('Reaction (r)')
#py.title('Cellulose')
#py.grid()
#py.show()
#
#py.figure(2)
#py.plot(t, CELL, label='CELL')
#py.legend(loc='best', numpoints=1)
#py.xlabel('Time (s)')
#py.ylabel('Species')
#py.title('Cellulose')
#py.grid()
#py.show()
