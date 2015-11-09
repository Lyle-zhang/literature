# Calculate terminal velocity, ut, from ut* and dp* for different particle sphericities
# from Kunii1991 book, pg 80-83, Eqs 31-33
# applicable for sphericity from 0.5 to 1.0

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

#--- Parameters

# air properties at T = 300K, P = 1 atm
rhog = 1.17     # density (kg/m^3)
ug = 1.85e-5    # dynamic viscosity (kg/ms)
g = 9.81        # gravity (m/s^2)

# particle properties
dp = 0.000207   # diameter of particle (m)
rhos = 2500     # density of particle (kg/m^3)
sp = 0.8        # sphericity of the particle, perfect sphere = 1.0

#--- Calculations

# dimensionless particle diameter, Kunii1991 Eq 31
dps = dp*(( (rhog*(rhos-rhog)*g) / (ug**2) )**(1/3))

# dimensionless terminal velocity, Kunii1991 Eq 33
uts = ( 18/(dps**2) + (2.335-1.744*sp)/(dps**0.5) )**-1

# terminal velocity, ut (m/s), Kunii1991 Eq 32
ut = uts*( (ug*(rhos-rhog)*g) / (rhog**2) )**(1/3)

# Reynolds number for particle, Re (-)
Re = rhog*ut*dp/ug

print('--- Kunii1991 ut*, dp* ---')
print('ut =', ut)
print('Re =', Re)
