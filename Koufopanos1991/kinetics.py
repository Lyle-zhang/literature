"""
Kinetic Reaction Scheme Functions for Fast Pyrolysis of Biomass.

Each function is for a particular kinetic scheme.
Reference for each scheme is provided as main author and publication year.
"""

# modules
# -----------------------------------------------------------------------------
import numpy as np

# Kinetics Function
# -----------------------------------------------------------------------------
    
def kn(T, B, C1, C2, rhow, dt, i, H):
    
    R = 0.008314 # universal gas constant, kJ/mol*K
    
    # A as pre-factor (1/s) and E as activation energy (kJ/mol)
    A1 = 9.973e-5;  G1 = 17254.4;   L1 = -9061227;  # biomass -> volatiles + gases
    A2 = 1.068e-3;  G2 = 10224.4;   L2 = -6123081;  # biomass -> char
    A3 = 5.7e5;     E3 = 81;        S = 1.45        # (vol+gases)1 -> (vol+gases)2
    
    # evaluate reaction rate constant for each reaction, 1/s
    K1 = A1 * np.exp((G1 / T[i]) + (L1 / T[i]**2))  # biomass -> volatiles + gases
    K2 = A2 * np.exp((G2 / T[i]) + (L2 / T[i]**2))  # biomass -> char
    K3 = A3 * np.exp(-E3 / (R * T[i]))              # (vol+gases)1 -> (vol+gases)2
    
    # reaction rates as mass fraction basis
    rB = -(K1+K2) * B[i-1]          # biomass rate
    rC1 = K2*B[i-1] - K3*C1[i-1]    # primary char rate
    rC2 = S*K3*C1[i-1]                # secondary char rate
    
    # update biomass and char mass fractions, (-)
    Bnew = B[i-1] + rB*dt
    C1new = C1[i-1] + rC1*dt
    C2new = C2[i-1] + rC2*dt
    
    # calculate heat of generation term
    rp = rhow*(rB+rC1+rC2)  # rate of pyrolysis
    g = H*rp                # heat generation
    
    # return the biomass and char mass fractions and heat of generation
    return Bnew, C1new, C2new, g
    
    