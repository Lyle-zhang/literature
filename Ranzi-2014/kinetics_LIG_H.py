# Lignin-H kinetics from Ranzi2013 paper
# plots species as function of time

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

#---- kinetic scheme from Table 1 Supplemental Material: Lignin (LIG-H)

# Arrhenius equation K = A exp(-E/RT) where
# A = pre-exponential factor, 1/s
# E = activation energy, kcal/kmol
# K = kinetic rate constant, 1/s

A11 = 0.67e13;       E11 = 37500
A14 = 0.5e8;         E14 = 30000
A15 = 33;            E15 = 15000
A16 = 2.4;           E16 = 12000
A17 = 0.4e9;         E17 = 30000
A18 = 0.083;         E18 = 8000

K11 = A11 * np.exp(-E11 / (R * T))   
K14 = A14 * np.exp(-E14 / (R * T))
K15 = A15 * np.exp(-E15 / (R * T))
K16 = A16 * T * np.exp(-E16 / (R * T))
K17 = A17 * np.exp(-E17 / (R * T))
K18 = A18 * T * np.exp(-E18 / (R * T))

# reactions in Table 1 Supplemental Material: Lignin (LIG-H)
# listed from top to bottom as 11, 14-18
# 11) LIG-H -> G1
# 14) LIG-OH -> G2
# 15) LIG-OH -> G3
# 16) LIG -> FE2MACR
# 17) LIG -> G4
# 18) LIG -> G5

# species array, sp, where rows = species, columns = time step
# sp[0] = LIG-H
# sp[1] = LIG-OH
# sp[2] = C3H6O
# sp[3] = LIG
# sp[4] = G{H2}
# sp[5] = H2O
# sp[6] = CH4
# sp[7] = CH3OH
# sp[8] = G{CH3OH}
# sp[9] = CO2
# sp[10] = CO
# sp[11] = G{CO}
# sp[12] = HCOOH
# sp[13] = G{COH2}
# sp[14] = G{CH4}
# sp[15] = G{C2H4}
# sp[16] = Char
# sp[17] = CH2O
# sp[18] = C2H4O

sp = np.zeros((19, nt))
sp[0, 0] = lig

for i in range(1, nt):
    sp[0, i] = sp[0, i-1] - K11*sp[0, i-1]*dt                                   # LIG-H
    sp[1, i] = sp[1, i-1] + K11*sp[0, i-1]*dt - (K14+K15)*sp[1, i-1]*dt         # LIG-OH
    sp[2, i] = sp[2, i-1] + K11*sp[0, i-1]*dt + K17*sp[3, i-1]*dt*0.2           # C3H6O
    sp[3, i] = sp[3, i-1] + K14*sp[1, i-1]*dt - (K16+K17+K18)*sp[3, i-1]*dt     # LIG
    sp[4, i] = sp[4, i-1] + K14*sp[1, i-1]*dt*0.15 + K15*sp[1, i-1]*dt*0.5      # G{H2}
    sp[5, i] = sp[5, i-1] + K14*sp[1, i-1]*dt*0.9 + K15*sp[1, i-1]*dt*1.5 + \
    K17*sp[3, i-1]*dt*0.95 + K18*sp[3, i-1]*dt*0.6                              # H2O
    sp[6, i] = sp[6, i-1] + K14*sp[1, i-1]*dt*0.1 + K15*sp[1, i-1]*dt*0.1 + \
    K17*sp[3, i-1]*dt*0.2 + K18*sp[3, i-1]*dt*0.2                               # CH4
    sp[7, i] = sp[7, i-1] + K14*sp[1, i-1]*dt*0.5 + K17*sp[3, i-1]*dt*0.4       # CH3OH
    sp[8, i] = sp[8, i-1] + K14*sp[1, i-1]*dt*0.5 + K15*sp[1, i-1]*dt*0.5 + \
    K18*sp[3, i-1]*dt*0.4                                                       # G{CH3OH}
    sp[9, i] = sp[9, i-1] + K14*sp[1, i-1]*dt*0.05                              # CO2
    sp[10, i] = sp[10, i-1] + K14*sp[1, i-1]*dt*0.3 + K15*sp[1, i-1]*dt*0.5 + \
    K17*sp[3, i-1]*dt + K18*sp[3, i-1]*dt*0.4                                   # CO
    sp[11, i] = sp[11, i-1] + K14*sp[1, i-1]*dt + K15*sp[1, i-1]*dt*1.6 + \
    K17*sp[3, i-1]*dt*0.45 + K18*sp[3, i-1]*dt*0.2                              # G{CO}
    sp[12, i] = sp[12, i-1] + K14*sp[1, i-1]*dt*0.05 + K17*sp[3, i-1]*dt*0.05   # HCOOH
    sp[13, i] = sp[13, i-1] + K14*sp[1, i-1]*dt*0.6 + K15*sp[1, i-1]*dt*3.9 + \
    K17*sp[3, i-1]*dt*0.5 + K18*sp[3, i-1]*dt*2                                 # G{COH2}
    sp[14, i] = sp[14, i-1] + K14*sp[1, i-1]*dt*0.35 + K15*sp[1, i-1]*dt*1.65 + \
    K17*sp[3, i-1]*dt*0.4 + K18*sp[3, i-1]*dt*0.4                               # G{CH4}
    sp[15, i] = sp[15, i-1] + K14*sp[1, i-1]*dt*0.2 + K15*sp[1, i-1]*dt*0.3 + \
    K17*sp[3, i-1]*dt*0.65 + K18*sp[3, i-1]*dt*0.5                              # G{C2H4}
    sp[16, i] = sp[16, i-1] + K14*sp[1, i-1]*dt*4.15 + K15*sp[1, i-1]*dt*10.15 + \
    K17*sp[3, i-1]*dt*5.5 + K18*sp[3, i-1]*dt*6                                 # Char
    sp[17, i] = sp[17, i-1] + K17*sp[3, i-1]*dt*0.2 + K18*sp[3, i-1]*dt*0.4     # CH2O
    sp[18, i] = sp[18, i-1] + K17*sp[3, i-1]*dt*0.2                             # C2H4O

#---- plot results as fraction vs t

py.figure(1)
py.plot(t, sp[0, :], label='LIG-H')
py.plot(t, sp[1, :], label='LIG-OH')
py.plot(t, sp[2, :], '--', label='C3H6O')
py.plot(t, sp[3, :], label='LIG')
py.plot(t, sp[4, :], label='G{H2}')
py.plot(t, sp[5, :], label='H2O')
py.plot(t, sp[6, :], label='CH4')
py.plot(t, sp[7, :], label='CH3OH')
py.plot(t, sp[8, :], label='G{CH3OH}')
py.plot(t, sp[9, :], label='CO2')
py.plot(t, sp[10, :], label='CO')
py.plot(t, sp[11, :], label='G{CO}')
py.plot(t, sp[12, :], label='HCOOH')
py.plot(t, sp[13, :], label='G{COH2}')
py.plot(t, sp[14, :], label='G{CH4}')
py.plot(t, sp[15, :], label='G{C2H4}')
py.plot(t, sp[16, :], label='Char')
py.plot(t, sp[17, :], label='CH2O')
py.plot(t, sp[18, :], label='C2H4O')
py.legend(loc='best', numpoints=1, prop={'size':10})
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Lignin (LIG-H) at T=%sK' % T)
py.grid()
py.show()

#py.figure(2)
#py.plot(t, sp[21, :], label='COUMARYL all')
#py.plot(t, sp[22, :], label='PHENOL all')
#py.plot(t, sp[23, :], label='H2O all')
#py.plot(t, sp[24, :], label='CO all')
#py.plot(t, sp[25, :], label='G{COH2} all')
#py.plot(t, sp[26, :], label='G{CH4} all')
#py.plot(t, sp[27, :], label='Char all')
#py.legend(loc='best', numpoints=1)
#py.xlabel('time (s)')
#py.ylabel('fraction (-)')
#py.title('Lignin (LIG-C) at T=%sK' % T)
#py.grid()
#py.show()
