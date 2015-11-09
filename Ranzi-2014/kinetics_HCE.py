# Hemicellulose kinetics from Ranzi2013 paper
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

#---- kinetic scheme from Table 1 Supplemental Material: Hemicellulose

# Arrhenius equation K = A exp(-E/RT) where
# A = pre-exponential factor, 1/s
# E = activation energy, kcal/kmol
# K = kinetic rate constant, 1/s

A5 = 0.33e10;       E5 = 31000
A6 = 1e9;           E6 = 32000
A7 = 0.05;           E7 = 8000
A8 = 0.9;           E8 = 11000
A9 = 0.33e10;       E9 = 33000

K5 = A5 * np.exp(-E5 / (R * T))   
K6 = A6 * np.exp(-E6 / (R * T))
K7 = A7 * T * np.exp(-E7 / (R * T))
K8 = A8 * np.exp(-E8 / (R * T))
K9 = A9 * np.exp(-E9 / (R * T))

# reactions in Table 1 Supplemental Material: Hemicellulose
# listed from top to bottom as 5-9
# 5) HCE -> HCE1, HCE2
# 6) HCE1 -> H20, CO2, HCOOH, CO, CH2O, C2H5OH, CH3OH, C2H4, GH2, GCO2, GCOH2
#            GCH3OH, GCH4, Char
# 7) HCE1 -> H2O, CO2, HCOOH, CO, GCO, GCO2, GCOH2, GCH4, GC2H4, Char
# 8) HCE1 -> XYLAN
# 9) HCE2 -> H2O, CO, CO2, CH2O, C2H5OH, HAA, HCOOH, GCH4, GCH3OH, GC2H4, GCO2
#            GCOH2, Char

# species array, sp, where rows = species, columns = time step
# sp[0] = HCE
# sp[1] = HCE1
# sp[2] = HCE2
# sp[3] = G2
# sp[4] = H2O
# sp[5] = CO2
# sp[6] = HCOOH
# sp[7] = CO
# sp[8] = CH2O
# sp[9] = C2H5OH
# sp[10] = CH3OH
# sp[11] = C2H4
# sp[12] = GH2
# sp[13] = GCO2
# sp[14] = GCOH2
# sp[15] = GCH3OH
# sp[16] = GCH4
# sp[17] = Char
# sp[18] = G3
# sp[19] = H2O
# sp[20] = CO2
# sp[21] = HCOOH
# sp[22] = CO
# sp[23] = GCO
# sp[24] = GCO2
# sp[25] = GCOH2
# sp[26] = GCH4
# sp[27] = GC2H4
# sp[28] = Char
# sp[29] = XYLAN
# sp[30] = G4
# sp[31] = H2O
# sp[32] = CO
# sp[33] = CO2
# sp[34] = CH2O
# sp[35] = C2H5OH
# sp[36] = HAA
# sp[37] = HCOOH
# sp[38] = GCH4
# sp[39] = GCH3OH
# sp[40] = GC2H4
# sp[41] = GCO2
# sp[42] = GCOH2
# sp[43] = Char
# sp[44] = H2O all
# sp[45] = CO2 all
# sp[46] = HCOOH all
# sp[47] = CO all
# sp[48] = CH2O all
# sp[49] = C2H5OH all
# sp[50] = GCO2 all
# sp[51] = GCOH2 all
# sp[52] = GCH3OH all
# sp[53] = GCH4 all
# sp[54] = Char all
# sp[55] = GC2H4 all

sp = np.zeros((56, nt))
sp[0, 0] = hemi

for i in range(1, nt):
    sp[0, i] = sp[0, i-1] - K5*sp[0, i-1]*dt                                    # HCE
    sp[1, i] = sp[1, i-1] + K5*sp[0, i-1]*dt*0.4 - (K6+K7+K8)*sp[1, i-1]*dt     # HCE1
    sp[2, i] = sp[2, i-1] + K5*sp[0, i-1]*dt*0.6 - K9*sp[2, i-1]*dt             # HCE2
    sp[3, i] = sp[3, i-1] + K6*sp[1, i-1]*dt                                    # G2
    sp[4, i] = sp[4, i-1] + K6*sp[1, i-1]*dt*0.025                              # H2O
    sp[5, i] = sp[5, i-1] + K6*sp[1, i-1]*dt*0.5                                # CO2
    sp[6, i] = sp[6, i-1] + K6*sp[1, i-1]*dt*0.025                              # HCOOH
    sp[7, i] = sp[7, i-1] + K6*sp[1, i-1]*dt*0.5                                # CO
    sp[8, i] = sp[8, i-1] + K6*sp[1, i-1]*dt*0.8                                # CH2O
    sp[9, i] = sp[9, i-1] + K6*sp[1, i-1]*dt*0.125                              # C2H5OH
    sp[10, i] = sp[10, i-1] + K6*sp[1, i-1]*dt*0.1                              # CH3OH
    sp[11, i] = sp[11, i-1] + K6*sp[1, i-1]*dt*0.25                             # C2H4
    sp[12, i] = sp[12, i-1] + K6*sp[1, i-1]*dt*0.125                            # GH2
    sp[13, i] = sp[13, i-1] + K6*sp[1, i-1]*dt*0.275                            # GCO2
    sp[14, i] = sp[14, i-1] + K6*sp[1, i-1]*dt*0.4                              # GCOH2
    sp[15, i] = sp[15, i-1] + K6*sp[1, i-1]*dt*0.45                             # GCH3OH
    sp[16, i] = sp[16, i-1] + K6*sp[1, i-1]*dt*0.325                            # GCH4
    sp[17, i] = sp[17, i-1] + K6*sp[1, i-1]*dt*0.875                            # Char
    sp[18, i] = sp[18, i-1] + K7*sp[1, i-1]*dt                                  # G3
    sp[19, i] = sp[19, i-1] + K7*sp[1, i-1]*dt*0.25                             # H2O
    sp[20, i] = sp[20, i-1] + K7*sp[1, i-1]*dt*0.5                              # CO2
    sp[21, i] = sp[21, i-1] + K7*sp[1, i-1]*dt*0.05                             # HCOOH
    sp[22, i] = sp[22, i-1] + K7*sp[1, i-1]*dt*0.3                              # CO
    sp[23, i] = sp[23, i-1] + K7*sp[1, i-1]*dt*0.15                             # GCO
    sp[24, i] = sp[24, i-1] + K7*sp[1, i-1]*dt*0.25                             # GCO2
    sp[25, i] = sp[25, i-1] + K7*sp[1, i-1]*dt*1.7                              # GCOH2
    sp[26, i] = sp[26, i-1] + K7*sp[1, i-1]*dt*0.625                            # GCH4
    sp[27, i] = sp[27, i-1] + K7*sp[1, i-1]*dt*0.375                            # GC2H4
    sp[28, i] = sp[28, i-1] + K7*sp[1, i-1]*dt*0.675                            # Char
    sp[29, i] = sp[29, i-1] + K8*sp[1, i-1]*dt                                  # XYLAN
    sp[30, i] = sp[30, i-1] + K9*sp[2, i-1]*dt                                  # G4
    sp[31, i] = sp[31, i-1] + K9*sp[2, i-1]*dt*0.2                              # H2O
    sp[32, i] = sp[32, i-1] + K9*sp[2, i-1]*dt*0.175                            # CO
    sp[33, i] = sp[33, i-1] + K9*sp[2, i-1]*dt*0.275                            # CO2
    sp[34, i] = sp[34, i-1] + K9*sp[2, i-1]*dt*0.5                              # CH2O
    sp[35, i] = sp[35, i-1] + K9*sp[2, i-1]*dt*0.1                              # C2H5OH
    sp[36, i] = sp[36, i-1] + K9*sp[2, i-1]*dt*0.2                              # HAA
    sp[37, i] = sp[37, i-1] + K9*sp[2, i-1]*dt*0.025                            # HCOOH
    sp[38, i] = sp[38, i-1] + K9*sp[2, i-1]*dt*0.25                             # GCH4
    sp[39, i] = sp[39, i-1] + K9*sp[2, i-1]*dt*0.3                              # GCH3OH
    sp[40, i] = sp[40, i-1] + K9*sp[2, i-1]*dt*0.275                            # GC2H4
    sp[41, i] = sp[41, i-1] + K9*sp[2, i-1]*dt*0.4                              # GCO2
    sp[42, i] = sp[42, i-1] + K9*sp[2, i-1]*dt*0.925                            # GCOH2
    sp[43, i] = sp[43, i-1] + K9*sp[2, i-1]*dt                                  # Char
    sp[44, i] = sp[4, i] + sp[19, i] + sp[31, i]                                # H2O all
    sp[45, i] = sp[5, i] + sp[20, i] + sp[33, i]                                # CO2 all
    sp[46, i] = sp[6, i] + sp[21, i] + sp[37, i]                                # HCOOH all
    sp[47, i] = sp[7, i] + sp[22, i] + sp[32, i]                                # CO all
    sp[48, i] = sp[8, i] + sp[34, i]                                            # CH2O all
    sp[49, i] = sp[9, i] + sp[35, i]                                            # C2H5OH all
    sp[50, i] = sp[13, i] + sp[24, i] + sp[41, i]                               # GCO2 all
    sp[51, i] = sp[14, i] + sp[25, i] + sp[42, i]                               # GCOH2 all
    sp[52, i] = sp[15, i] + sp[39, i]                                           # GCH3OH all
    sp[53, i] = sp[16, i] + sp[26, i] + sp[38, i]                               # GCH4 all
    sp[54, i] = sp[17, i] + sp[28, i] + sp[43, i]                               # Char all
    sp[55, i] = sp[27, i] + sp[40, i]                                           # GC2H4 all
    

#---- plot results as fraction vs t

py.figure(1)
py.plot(t, sp[0, :], label='HCE')
py.plot(t, sp[1, :], label='HCE1')
py.plot(t, sp[2, :], label='HCE2')
py.plot(t, sp[3, :], label='G2')
py.plot(t, sp[18, :], label='G3')
py.plot(t, sp[30, :], label='G4')
py.legend(loc='best', numpoints=1)
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Hemicellulose at T=%sK' % T)
py.grid()
py.show()

py.figure(2)
py.plot(t, sp[4, :], label='H2O')
py.plot(t, sp[5, :], label='CO2')
py.plot(t, sp[6, :], '--', label='HCOOH')
py.plot(t, sp[7, :], '--', label='CO')
py.plot(t, sp[8, :], label='CH2O')
py.plot(t, sp[9, :], label='C2H5OH')
py.plot(t, sp[10, :], label='CH3OH')
py.plot(t, sp[11, :], label='C2H4')
py.plot(t, sp[12, :], '--', label='GH2')
py.plot(t, sp[13, :], label='GCO2')
py.plot(t, sp[14, :], label='GCOH2')
py.plot(t, sp[15, :], label='GCH3OH')
py.plot(t, sp[16, :], label='GCH4')
py.plot(t, sp[17, :], label='Char')
py.legend(loc='best', numpoints=1, prop={'size':12})
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Hemicellulose (G2) at T=%sK' % T)
py.grid()
py.show()

py.figure(3)
py.plot(t, sp[19, :], label='H2O')
py.plot(t, sp[20, :], label='CO2')
py.plot(t, sp[21, :], label='HCOOH')
py.plot(t, sp[22, :], label='CO')
py.plot(t, sp[23, :], label='GCO')
py.plot(t, sp[24, :], '--', label='GCO2')
py.plot(t, sp[25, :], label='GCOH2')
py.plot(t, sp[26, :], label='GCH4')
py.plot(t, sp[27, :], label='GC2H4')
py.plot(t, sp[28, :], label='Char')
py.plot(t, sp[29, :], label='XYLAN')
py.legend(loc='best', numpoints=1, prop={'size':12})
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Hemicellulose (G3) at T=%sK' % T)
py.grid()
py.show()

py.figure(4)
py.plot(t, sp[31, :], label='H2O')
py.plot(t, sp[32, :], label='CO')
py.plot(t, sp[33, :], label='CO2')
py.plot(t, sp[34, :], label='CH2O')
py.plot(t, sp[35, :], label='C2H5OH')
py.plot(t, sp[36, :], '--', label='HAA')
py.plot(t, sp[37, :], label='HCOOH')
py.plot(t, sp[38, :], label='GCH4')
py.plot(t, sp[39, :], label='GCH3OH')
py.plot(t, sp[40, :], label='GC2H4')
py.plot(t, sp[41, :], label='GCO2')
py.plot(t, sp[42, :], label='GCOH2')
py.plot(t, sp[43, :], label='Char')
py.legend(loc='best', numpoints=1, prop={'size':12})
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Hemicellulose (G4) at T=%sK' % T)
py.grid()
py.show()

py.figure(5)
py.plot(t, sp[44, :], label='H2O')
py.plot(t, sp[45, :], label='CO2')
py.plot(t, sp[46, :], label='HCOOH')
py.plot(t, sp[47, :], label='CO')
py.plot(t, sp[48, :], label='CH2O')
py.plot(t, sp[49, :], label='C2H5OH')
py.plot(t, sp[50, :], label='GCO2')
py.plot(t, sp[51, :], label='GCOH2')
py.plot(t, sp[52, :], label='GCH3OH')
py.plot(t, sp[53, :], label='GCH4')
py.plot(t, sp[54, :], label='Char')
py.plot(t, sp[55, :], label='GC2H4')
py.legend(loc='best', numpoints=1, prop={'size':12})
py.xlabel('time (s)')
py.ylabel('fraction (-)')
py.title('Hemicellulose (all) at T=%sK' % T)
py.grid()
py.show()