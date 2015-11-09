# Compare equations in Santos book to the equations that he referenced
# Created: 2/6/2014
# Author: Gavin Wiggins

# use Python 3 print function
from __future__ import print_function

# libraries and packages
import numpy as np

# -----------------------------------------------------------------------------
# example parameters

umf = 0.063     # minimum fluidization velocity, m/s
u = 3*umf       # gas superficial velocity, m/s
dd = 0.055      # bed diameter, m
zst = 0.085     # static bed height, m
g = 9.81        # gravity, m/s^2
s = (np.pi*dd**2)/4.0  # bed cross-sectional area, m^2

# Eq 14.28b - minimum bubble diameter for porous plate for diameter at z = 0
dbmin = 3.77*(u-umf)**2/9.81

# Eq 14.27 - maximum bubble diameter
dbmax = 2.59*(g**-0.2)*(s*(u-umf))**0.4

# Eq 14.26 - bubble diameter according to vertical position z
db = dbmax - (dbmax - dbmin)*np.exp(-0.3*(zst/dd))
    
# -----------------------------------------------------------------------------
# Mori and Ken equations for Jack Halow case
# *note that these are in cm not meters

dbed = 5.5  				  # bed diameter, cm
areabed = (np.pi*dbed**2)/4.0   # bed area, cm^2
umfM = 6.3                      # minimum fluidization velocity, cm/s
uM = 3*umfM                     # gas superficial velocity, cm/s
h = 8.5                         # height of bed, cm

# max bubble diameter, cm
dbmaxM = 0.652*(areabed*(uM-umfM))**(2/5.0)

# initial bubble diameter for porous plate, cm
dbminM = 0.00376*(uM-umfM)**2

# bubble growth correlation for bubble diameter, cm
dbM = dbmaxM - (dbmaxM - dbminM)*np.exp(-0.3*(h/dbed))

# -----------------------------------------------------------------------------
# Print results - note that Mori must be converted to m from cm
print('\nSantos dbmin =',dbmin)
print('Mori dbmin =',dbminM/100)

print('\nSantos dbmax =',dbmax)
print('Mori dbmax =',dbmaxM/100)

print('\nSantos db(0.085) =',db)
print('Mori db(0.085) =',dbM/100)

# looks like the results match the referenced equations' results