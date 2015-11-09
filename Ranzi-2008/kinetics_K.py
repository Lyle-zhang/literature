# Plots reaction rate constant, K, as function of temperature, T from
# cellulose, hemicellulose, and lignin components.
# kinetics from Ranzi2008 paper

# use Python 3 print function
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

#---- global parameters

R = 1.987   # gas constant, kcal/kmol*K

T = np.linspace(700, 900)  # temperature range, K

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

#---- plot T vs K for cellulose, hemicellulose, lignin

py.figure(1)
py.plot(T, K1, label='K1')
py.plot(T, K2, label='K2')
py.plot(T, K3, label='K3')
py.plot(T, K4, label='K4')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate (1/s)')
py.title('Cellulose')
#py.yscale('log')
py.grid()
py.show()

py.figure(2)
py.plot(T, K5, label='K5')
py.plot(T, K6, label='K6')
py.plot(T, K7, label='K7')
py.plot(T, K8, label='K8')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate (1/s)')
py.title('Hemicellulose')
#py.yscale('log')
py.grid()
py.show()

py.figure(3)
py.plot(T, K9, label='K9')
py.plot(T, K10, label='K10')
py.plot(T, K11, label='K11')
py.plot(T, K12, label='K12')
py.plot(T, K13, label='K13')
py.plot(T, K14, label='K14')
py.plot(T, K15, label='K15')
py.plot(T, K16, label='K16')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate (1/s)')
py.title('Lignin')
#py.yscale('log')
py.grid()
py.show()
