# Horio and Nonaka bubble diameter correlation
# see Santos2010 book, Eq 14.30, pg 322

# use Python 3 print function
from __future__ import print_function

# libraries and packages
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as py

#------------------------------------------------------------------------------
# parameters

umf = 0.063         # minimum fluidization velocity, m/s
dbed = 0.055        # bed diameter, m
z0 = 0              # position bubbles are generated, m
z = 0.046           # bed vertical position, m

g = 9.81            # gravity, m/s^2

#------------------------------------------------------------------------------
# calculations

m = 3                       # multiplier for Umf
u = m*umf                   # gas superficial velocity, m/s

abed = (np.pi*dbed**2)/4.0  # bed cross-sectional area, m^2

# Horio & Nonaka correlations
# general form of equation: (term1)^power1 * (term2)^power2 = term3

dbmax = 2.59*(g**-0.2)*(abed*(u-umf))**0.4
dbmin = 3.77*(u-umf)**2/g

c1 = 2.56*10**-2*((dbed / g)**0.5/umf)

c2 = (c1**2 + (4*dbmax)/dbed)**0.5

c3 = 0.25*dbed*(c1 + c2)**2

dbeq = 0.25*dbed*(-c1 + (c1**2 + 4*(dbmax/dbed))**0.5 )**2

power1 = 1 - c1/c2

power2 = 1 + c1/c2

term3 = np.exp(-0.3*(z - z0)/dbed)

def dB(d):
    term1 = (np.sqrt(d) - np.sqrt(dbeq)) / (np.sqrt(dbmin) - np.sqrt(dbeq))
    term2 = (np.sqrt(d) + np.sqrt(c3)) / (np.sqrt(dbmin) + np.sqrt(c3))
    return term1**power1 * term2**power2 - term3

dbub = opt.newton_krylov(dB, [1e-6, 0.055])
print('dbub = ', dbub)

# solve bubble diameter for range of bed height z

db = []
zz = []

for i in np.arange(0, 0.118, 0.001):
    z = i
    zz.append(i)
    term3 = np.exp(-0.3*(z - z0)/dbed)
    dbub = opt.newton_krylov(dB, [1e-6, 0.055])
    db.append(dbub[0])

py.figure()
py.plot(zz, db, '-g')
py.title('Bubble Characteristics: ' + str(m) + '*Umf')
py.grid()
py.rcParams['xtick.major.pad'] = 6
py.rcParams['ytick.major.pad'] = 6
py.show()
