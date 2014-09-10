# Lignin-C kinetics from Ranzi2013 paper
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

#---- kinetic scheme from Table 1 Supplemental Material: Lignin (LIG-C)

# Arrhenius equation K = A exp(-E/RT) where
# A = pre-exponential factor, 1/s
# E = activation energy, kcal/kmol
# K = kinetic rate constant, 1/s

A10 = 1.33e15;       E10 = 48500
A13 = 1.6e6;         E13 = 31500

K10 = A10 * np.exp(-E10 / (R * T))   
K13 = A13 * np.exp(-E13 / (R * T))

# reactions in Table 1 Supplemental Material: Lignin (LIG-C)
# listed from top to bottom as 10, 13
# 10) LIG-C -> G1
# 13) LIG-CC -> G2

# species array, sp, where rows = species, columns = time step
# sp[0] = LIG-C
# sp[1] = LIG-CC
# sp[2] = COUMARYL
# sp[3] = PHENOL
# sp[4] = C2H4
# sp[5] = H2O
# sp[6] = CH2O
# sp[7] = CO
# sp[8] = G{COH2}
# sp[9] = G{CH4}
# sp[10] = Char
# sp[11] = COUMARYL
# sp[12] = PHENOL
# sp[13] = HAA
# sp[14] = H2O
# sp[15] = CO
# sp[16] = G{CH4}
# sp[17] = G{C2H4}
# sp[18] = G{COH2}
# sp[19] = G{CO}
# sp[20] = Char
# sp[21] = COUMARYL all
# sp[22] = PHENOL all
# sp[23] = H2O all
# sp[24] = CO all
# sp[25] = G{COH2} all
# sp[26] = G{CH4} all
# sp[27] = Char all

sp = np.zeros((28, nt))
sp[0, 0] = lig

for i in range(1, nt):
    sp[0, i] = sp[0, i-1] - K10*sp[0, i-1]*dt                           # LIG-C
    sp[1, i] = sp[1, i-1] + K10*sp[0, i-1]*dt*0.35 - K13*sp[1, i-1]*dt  # LIG-CC
    sp[2, i] = sp[2, i-1] + K10*sp[0, i-1]*dt*0.1                       # COUMARYL
    sp[3, i] = sp[3, i-1] + K10*sp[0, i-1]*dt*0.08                      # PHENOL
    sp[4, i] = sp[4, i-1] + K10*sp[0, i-1]*dt*0.41                      # C2H4
    sp[5, i] = sp[5, i-1] + K10*sp[0, i-1]*dt                           # H2O
    sp[6, i] = sp[6, i-1] + K10*sp[0, i-1]*dt*0.3                       # CH2O
    sp[7, i] = sp[7, i-1] + K10*sp[0, i-1]*dt*0.32                      # CO
    sp[8, i] = sp[8, i-1] + K10*sp[0, i-1]*dt*0.7                       # G{COH2}
    sp[9, i] = sp[9, i-1] + K10*sp[0, i-1]*dt*0.495                     # G{CH4}
    sp[10, i] = sp[10, i-1] + K10*sp[0, i-1]*dt*5.735                   # Char
    sp[11, i] = sp[11, i-1] + K13*sp[1, i-1]*dt*0.3                     # COUMARYL
    sp[12, i] = sp[12, i-1] + K13*sp[1, i-1]*dt*0.2                     # PHENOL
    sp[13, i] = sp[13, i-1] + K13*sp[1, i-1]*dt*0.35                    # HAA
    sp[14, i] = sp[14, i-1] + K13*sp[1, i-1]*dt*0.7                     # H2O
    sp[15, i] = sp[15, i-1] + K13*sp[1, i-1]*dt*0.4                     # CO
    sp[16, i] = sp[16, i-1] + K13*sp[1, i-1]*dt*0.65                    # G{CH4}
    sp[17, i] = sp[17, i-1] + K13*sp[1, i-1]*dt*0.6                     # G{C2H4}
    sp[18, i] = sp[18, i-1] + K13*sp[1, i-1]*dt                         # G{COH2}
    sp[19, i] = sp[19, i-1] + K13*sp[1, i-1]*dt*0.4                     # G{CO}
    sp[20, i] = sp[20, i-1] + K13*sp[1, i-1]*dt*6.75                    # Char
    sp[21, i] = sp[21, i-1] + K10*sp[0, i-1]*dt*0.1 + K13*sp[1, i-1]*dt*0.3     # COUMARYL all
    #sp[22, i] = sp[2, i] + sp[11, i]                                           # COUMARYL all, another approach
    sp[22, i] = sp[22, i-1] + K10*sp[0, i-1]*dt*0.08 + K13*sp[1, i-1]*dt*0.2    # PHENOL all
    sp[23, i] = sp[23, i-1] + K10*sp[0, i-1]*dt + K13*sp[1, i-1]*dt*0.7         # H2O all
    sp[24, i] = sp[24, i-1] + K10*sp[0, i-1]*dt*0.32 + K13*sp[1, i-1]*dt*0.4    # CO all
    sp[25, i] = sp[25, i-1] + K10*sp[0, i-1]*dt*0.7 + K13*sp[1, i-1]*dt         # G{COH2} all
    sp[26, i] = sp[26, i-1] + K10*sp[0, i-1]*dt*0.495 + K13*sp[1, i-1]*dt*0.65  # G{CH4} all
    sp[27, i] = sp[27, i-1] + K10*sp[0, i-1]*dt*5.735 + K13*sp[1, i-1]*dt*6.75  # Char all

#---- plot results as fraction vs t

py.figure(1)
py.plot(t, sp[0, :], label='LIG-C')
py.plot(t, sp[1, :], label='LIG-CC')
py.plot(t, sp[2, :], label='COUMARYL')
py.plot(t, sp[3, :], label='PHENOL')
py.plot(t, sp[4, :], label='C2H4')
py.plot(t, sp[5, :], label='H2O')
py.plot(t, sp[6, :], label='CH2O')
py.plot(t, sp[7, :], label='CO')
py.plot(t, sp[8, :], label='G{COH2}')
py.plot(t, sp[9, :], label='G{CH4}')
py.plot(t, sp[10, :], label='Char')
py.plot(t, sp[11, :], label='COUMARYL')
py.plot(t, sp[12, :], label='PHENOL')
py.plot(t, sp[13, :], label='HAA')
py.plot(t, sp[14, :], label='H2O')
py.plot(t, sp[15, :], label='CO')
py.plot(t, sp[16, :], label='G{CH4}')
py.plot(t, sp[17, :], label='G{C2H4}')
py.plot(t, sp[18, :], label='G{COH2}')
py.plot(t, sp[19, :], label='G{CO}')
py.plot(t, sp[20, :], label='Char')
py.legend(loc='best', numpoints=1, prop={'size':10})
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Lignin (LIG-C) at T=%sK' % T)
py.grid()
py.show()

py.figure(2)
py.plot(t, sp[21, :], label='COUMARYL all')
py.plot(t, sp[22, :], label='PHENOL all')
py.plot(t, sp[23, :], label='H2O all')
py.plot(t, sp[24, :], label='CO all')
py.plot(t, sp[25, :], label='G{COH2} all')
py.plot(t, sp[26, :], label='G{CH4} all')
py.plot(t, sp[27, :], label='Char all')
py.legend(loc='best', numpoints=1)
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Lignin (LIG-C) at T=%sK' % T)
py.grid()
py.show()
