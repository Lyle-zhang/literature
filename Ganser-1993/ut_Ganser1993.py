# Calculate terminal velocity, ut, from drag coefficient, Cd, and particle sphericity
# Cd from Ganser1993, Eq 18
# Cd also referenced in Cui2007, Eq 2
# Cd also referenced in Chhabra1999, Eq 6

# according to Chhabra1999, the Ganser1993 Cd correlation is applicable for
# sphericity values from 0.09 to 1.0

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries
import numpy as np

#---- Parameters

# air properties at T = 300K, P = 1 atm
rhog = 1.17     # density (kg/m^3)
ug = 1.85e-5    # dynamic viscosity (kg/ms)
g = 9.81        # gravity (m/s^2)

# particle properties
dp = 0.000207   # diameter of particle (m)
rhos = 2500     # density of particle (kg/m^3)
sp = 0.8        # sphericity of the particle, perfect sphere = 1.0

#---- Calculations

# papers Cui2007 and Chhabra1999 leave out the -2.25*dv/D term
K1 = (1/3 + 2/3*(sp**-0.5))**(-1)           # Stokes' shape factor
K2 = 10**(1.8148*((-np.log(sp))**0.5743))   # Newton's shape factor

# guess range for terminal velocity, ut (m/s)
ut = np.arange(0.001, 20, 0.001)           

# Re, Reynolds number (-)
# Cd, Ganser1993, drag coefficient as function of Re, K1, K2 (-)
# Cdd, drag coefficient as function of ut, etc. (-)
Re = (dp*rhog*ut)/ug
Cd = (24/(Re*K1))*(1 + 0.1118*((Re*K1*K2)**0.6567)) + (0.4305*K2)/(1+(3305/(Re*K1*K2)))
Cdd = (4*g*dp*(rhos-rhog))/(3*(ut**2)*rhog)

delta = np.abs(Cd-Cdd)  # compare difference between Cd and Cdd
idx = np.argmin(delta)  # find index of minimum value in delta

# print results to console
print('--- Ganser1993 ---')
print('ut =', ut[idx])
print('Re =', Re[idx])
print('Cd =', Cd[idx])
print('Cdd =', Cdd[idx])
print('delta =', delta[idx])