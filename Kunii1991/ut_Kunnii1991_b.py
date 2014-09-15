# Calculate terminal velocity, ut, from drag coefficient, Cd, and particle sphericity
# from Kunii1991 book, pg 80, Eqs 28-29

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries
import numpy as np

#--- INPUTS

# air properties at T = 300K, P = 1 atm
rhog = 1.17     # density (kg/m^3)
ug = 1.85e-5    # dynamic viscosity (kg/ms)
g = 9.81        # gravity (m/s^2)

# particle properties
dp = 0.000207   # diameter of particle (m)
rhos = 2500     # density of particle (kg/m^3)
sp = 0.8        # sphericity of the particle, perfect sphere = 1.0 

#--- OUTPUTS

# guess range for terminal velocity, ut (m/s)
ut = np.arange(0.001, 20, 0.001) 

# Re, Reynolds number (-)
# Cd, drag coefficient as function of Re and sphericity (-)
# Cdd, drag coefficient as function of ut, etc. (-)
Re = (dp*rhog*ut)/ug
Cd = (24/Re)*(1+(8.1716*np.exp(-4.0655*sp))*Re**(0.0964+0.5565*sp)) + (73.69*np.exp(-5.0748*sp)*Re)/(Re+5.378*np.exp(6.2122*sp))
Cdd = (4*g*dp*(rhos-rhog))/(3*(ut**2)*rhog)

delta = np.abs(Cd-Cdd)  # compare difference between Cd and Cdd
idx = np.argmin(delta)  # find index of minimum value in delta

# print results to console
print('--- Kunii1991 ---')
print('ut =', ut[idx])
print('Re =', Re[idx])
print('Cd =', Cd[idx])
print('Cdd =', Cdd[idx])
print('delta =', delta[idx])

