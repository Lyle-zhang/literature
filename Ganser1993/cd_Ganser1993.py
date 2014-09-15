# Drag coefficient, Cd, for nonspherical particles
# see pg 150 of Ganser1993 article, Eq 18 and Table 7
# by G.W. 4/2/2014

# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries
import numpy as np
import matplotlib.pyplot as py

#--- parameters for Jack Halow case
dp = 0.000207   # bed particle diameter (m)
ps = 2500       # density of bed particle (kg/m^3)
pg = 1.17       # density of air (kg/m^3) at T = 300K, P = 1atm
ug = 1.85e-5    # dynamic viscosity of air (kg/ms) at T = 300K, P = 1atm
g = 9.81        # gravity (m/s^2)
sp = 0.8        # sphericity of the particle, perfect sphere = 1.0

#--- calculate drag coefficient, Cd

U = 1.89    # velocity (m/s)

Re = (dp*pg*U) / ug     # Reynolds number

K1 = (1/3 + 2/3*(sp**-0.5))**(-1)   # Stokes' shape factor

K2 = 10**(1.8148*((-np.log(sp))**0.5743))   # Newton's shape factor

# drag coefficient
Cd = (24*K2)/(Re*K1*K2)*(1+0.1118*(Re*K1*K2)**0.6567) + (0.4305*K2)/(1+(3305/(Re*K1*K2)))

print('Cd =', Cd)

#--- calculate Cd from range of Re
# note that the only inputs needed are Re and sphericity (K1, K2)
Ree = []
Cdd = []

i = np.logspace(-1, 4)
Ree.append(i)
Cdd.append((24*K2)/(i*K1*K2)*(1+0.1118*(i*K1*K2)**0.6567) + (0.4305*K2)/(1+(3305/(i*K1*K2))))

# plot Re vs Cd
# this plot should be similar to Fig 2 in Ganser1993 paper

py.figure(1)
py.loglog(Ree[0], Cdd[0])
py.title('Sphericity = %.1f' % sp)
py.xlabel('Reynolds number, Re')
py.ylabel('Drag coefficient, Cd')
py.grid(True, which='both', ls='-', alpha=0.4)
py.show()