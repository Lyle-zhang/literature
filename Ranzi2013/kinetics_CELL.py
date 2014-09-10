# Cellulose kinetics from Ranzi2013 paper
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
T = 773         # temperature, K

dt = 0.001                      # time step, s
tmax = 1                        # time max, s
t = np.arange(0, tmax+dt, dt)   # time range, s
nt = len(t)                     # number of time steps

#---- kinetic scheme from Table 1 Supplemental Material: Cellulose

# A = pre-exponential factor, 1/s and E = activation energy, kcal/kmol
A1 = 4e13;      E1 = 45000
A2 = 0.5e9;     E2 = 29000
A3 = 1.8;       E3 = 10000
A4 = 4e7;       E4 = 31000

# kinetic rate contants, Arrhenius equation K = A exp(-E/RT)
K1 = A1 * np.exp(-E1 / (R * T))   
K2 = A2 * np.exp(-E2 / (R * T))
K3 = A3 * T * np.exp(-E3 / (R * T))
K4 = A4 * np.exp(-E4 / (R * T))

# species array, sp, where rows = species, columns = time step
# sp[0] = CELL          cellulose
# sp[1] = G1            group 1
# sp[2] = H2O_1         H2O in group 1
# sp[3] = Char_1        Char in group 1
# sp[4] = CELLA         cellulose active
# sp[5] = LVG           levoglucosan
# sp[6] = G2            group 2
# sp[7] = HAA
# sp[8] = GLYOX
# sp[9] = C2H4O
# sp[10] = HMFU
# sp[11] = C3H6O
# sp[12] = CO2
# sp[13] = H2
# sp[14] = CH2O
# sp[15] = CO
# sp[16] = CH4
# sp[17] = H2O_2        H2O in group 2
# sp[18] = HCOOH
# sp[19] = Char_2       Char in group 2
# sp[20] = H2O_all      H2O = H2O_1 + H2O_2
# sp[21] = Char_all     Char = Char_1 + Char_2

sp = np.zeros((22, nt))
sp[0, 0] = cell

for i in range(1, nt):
    sp[0, i] = sp[0, i-1] - (K1 + K4)*sp[0, i-1]*dt                     # CELL
    sp[1, i] = sp[1, i-1] + K4*sp[0, i-1]*dt                            # G1
    sp[2, i] = sp[1, i] * 5                                             # H2O_1
    sp[3, i] = sp[1, i] * 6                                             # Char_1
    sp[4, i] = sp[4, i-1] + K1*sp[0, i-1]*dt - (K2 + K3)*sp[4, i-1]*dt  # CELLA
    sp[5, i] = sp[5, i-1] + K3*sp[4, i-1]*dt                            # LVG
    sp[6, i] = sp[6, i-1] + K2*sp[4, i-1]*dt                            # G2
    sp[7, i] = sp[6, i] * 0.8                                           # HAA
    sp[8, i] = sp[6, i] * 0.2                                           # GLYOX
    sp[9, i] = sp[6, i] * 0.1                                           # C2H4O
    sp[10, i] = sp[6, i] * 0.25                                         # HMFU
    sp[11, i] = sp[6, i] * 0.3                                          # C3H6O
    sp[12, i] = sp[6, i] * 0.21                                         # CO2
    sp[13, i] = sp[6, i] * 0.1                                          # H2
    sp[14, i] = sp[6, i] * 0.4                                          # CH2O
    sp[15, i] = sp[6, i] * 0.16                                         # CO
    sp[16, i] = sp[6, i] * 0.1                                          # CH4
    sp[17, i] = sp[6, i] * 0.83                                         # H2O_2
    sp[18, i] = sp[6, i] * 0.02                                         # HCOOH
    sp[19, i] = sp[6, i] * 0.61                                         # Char_2
    sp[20, i] = sp[2, i] + sp[17, i]                                    # H2O_all
    sp[21, i] = sp[3, i] + sp[19, i]                                    # Char_all


#---- plot results as fraction vs t

py.figure(1)
py.plot(t, sp[0, :], label='CELL')
py.plot(t, sp[4, :], label='CELLA')
py.plot(t, sp[1, :], label='G1')
py.plot(t, sp[20, :], label='H2O_all')
py.plot(t, sp[21, :], label='Char_all')
py.plot(t, sp[6, :], label='G2')
py.legend(loc='best', numpoints=1)
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Cellulose at T=%sK' % T)
py.grid()
py.show()

py.figure(2)
py.plot(t, sp[2, :], label='H2O')
py.plot(t, sp[3, :], label='Char')
py.plot(t, sp[5, :], label='LVG')
py.plot(t, sp[7, :], label='HAA')
py.plot(t, sp[8, :], label='GLYOX')
py.plot(t, sp[9, :], label='C2H4O')
py.plot(t, sp[10, :], label='HMFU')
py.plot(t, sp[11, :], label='C3H6O')
py.legend(loc='best', numpoints=1)
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Cellulose at T=%sK' % T)
py.grid()
py.show()

py.figure(3)
py.plot(t, sp[12, :], label='CO2')
py.plot(t, sp[13, :], label='H2')
py.plot(t, sp[14, :], label='CH2O')
py.plot(t, sp[15, :], label='CO')
py.plot(t, sp[16, :], label='CH4')
py.plot(t, sp[17, :], label='H2O')
py.plot(t, sp[18, :], label='HCOOH')
py.plot(t, sp[19, :], label='Char')
py.legend(loc='best', numpoints=1)
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Cellulose at T=%sK' % T)
py.grid()
py.show()

py.figure(4)
py.plot(t, sp[2, :], 'b-.', label='H2O_1')
py.plot(t, sp[3, :], 'r-.', label='Char_1')
py.plot(t, sp[17, :], 'b--', label='H2O_2')
py.plot(t, sp[19, :], 'r--', label='Char_2')
py.plot(t, sp[20, :], 'b', label='H2O_all')
py.plot(t, sp[21, :], 'r', label='Char_all')
py.legend(loc='best', numpoints=1)
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Cellulose at T=%sK' % T)
py.grid()
py.show()