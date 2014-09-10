# Kinetics from Ranzi2013 paper
# plots reaction rate constant, K, as function of temperature, T

# use Python 3 print function
from __future__ import print_function
from __future__ import division

import numpy as np
import matplotlib.pyplot as py

#---- global parameters

R = 1.987                   # gas constant, kcal/kmol*K
T = np.linspace(700, 900)   # temperature range, K

#---- Cellulose: kinetic scheme from Table 1 Supplemental Material

# A = pre-exponential factor, 1/s
# E = activation energy, kcal/kmol
A1 = 4e13;      E1 = 45000
A2 = 0.5e9;     E2 = 29000
A3 = 1.8;       E3 = 10000
A4 = 4e7;       E4 = 31000

# reaction rate constants, Arrhenius equation K = A*exp(-E/RT)
K1 = A1 * np.exp(-E1 / (R * T))   
K2 = A2 * np.exp(-E2 / (R * T))
K3 = A3 * T * np.exp(-E3 / (R * T))
K4 = A4 * np.exp(-E4 / (R * T))

#---- Hemicellulose: kinetic scheme from Table 1 Supplemental Material

A5 = 0.33e10;       E5 = 31000
A6 = 1e9;           E6 = 32000
A7 = 0.05;          E7 = 8000
A8 = 0.9;           E8 = 11000
A9 = 0.33e10;       E9 = 33000

K5 = A5 * np.exp(-E5 / (R * T))   
K6 = A6 * np.exp(-E6 / (R * T))
K7 = A7 * T * np.exp(-E7 / (R * T))
K8 = A8 * np.exp(-E8 / (R * T))
K9 = A9 * np.exp(-E9 / (R * T))

#---- Lignin: kinetic scheme from Table 1 Supplemental Material

A10 = 1.33e15;       E10 = 48500
A11 = 0.67e13;       E11 = 37500
A12 = 0.33e9;        E12 = 25500 
A13 = 1.6e6;         E13 = 31500
A14 = 0.5e8;         E14 = 30000
A15 = 33;            E15 = 15000
A16 = 2.4;           E16 = 12000
A17 = 0.4e9;         E17 = 30000
A18 = 0.083;         E18 = 8000

K10 = A10 * np.exp(-E10 / (R * T))   
K11 = A11 * np.exp(-E11 / (R * T))
K12 = A12 * np.exp(-E12 / (R * T))
K13 = A13 * np.exp(-E13 / (R * T))      
K14 = A14 * np.exp(-E14 / (R * T))   
K15 = A15 * np.exp(-E15 / (R * T))   
K16 = A16 * T * np.exp(-E16 / (R * T))   
K17 = A17 * np.exp(-E17 / (R * T))   
K18 = A18 * T * np.exp(-E18 / (R * T))

#---- plot T vs K for cellulose, hemicellulose, lignin

py.figure(1)
py.plot(T, K1, label=r'K1, CELL $\to$ CELLA')
py.plot(T, K2, label=r'K2, CELLA $\to$ G2')
py.plot(T, K3, label=r'K3, CELLA $\to$ LVG')
py.plot(T, K4, label=r'K4, CELL $\to$ G1')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate (1/s)')
py.title('Cellulose')
#py.yscale('log')
py.grid()
py.show()

py.figure(2)
py.plot(T, K5, label=r'K5, HCE $\to$ G1')
py.plot(T, K6, label=r'K6, HCE1 $\to$ G2')
py.plot(T, K7, label=r'K7, HCE1 $\to$ G3')
py.plot(T, K8, label=r'K8, HCE1 $\to$ Xylan')
py.plot(T, K9, label=r'K9, HCE2 $\to$ G4')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate (1/s)')
py.title('Hemicellulose')
#py.yscale('log')
py.grid()
py.show()

py.figure(3)
py.plot(T, K10, label=r'K10, LIG-C $\to$ G1')
py.plot(T, K11, label=r'K11, LIG-H $\to$ G1')
py.plot(T, K12, label=r'K12, LIG-O $\to$ G1')
py.plot(T, K13, label=r'K13, LIG-CC $\to$ G2')
py.plot(T, K14, label=r'K14, LIG-OH $\to$ G2')
py.plot(T, K15, label=r'K15, LIG-OH $\to$ G3')
py.plot(T, K16, label=r'K16, LIG $\to$ FE2MACR')
py.plot(T, K17, label=r'K17, LIG $\to$ G4')
py.plot(T, K18, label=r'K18, LIG $\to$ G5')
py.legend(loc='best', numpoints=1)
py.xlabel('Temperature (K)')
py.ylabel('Reaction Rate (1/s)')
py.title('Lignin')
#py.yscale('log')
py.grid()
py.show()
