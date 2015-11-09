# Fluidized bed calculations for fbexp, zexp

# applied to Jack Halow bubbling fluidized bed case
# equations from Santos2010 [citations from book in brackets]

# Created: 2/19/2014
# Author: Gavin Wiggins
#==============================================================================

# use Python 3 print function
from __future__ import print_function

# fluidized bed parameters

umf = 0.063                 # minimum fluidization velocity, m/s
u = 3.0*umf                 # gas superficial velocity, m/s
dbed = 0.055                # bed diameter, m
zmf = 0.085                 # bed height at minimum fluidizaiton, m
emf = 0.48                  # void fraction at minimum fluidization
rhop = 2500                 # particle density, kg/m^3

g = 9.81                    # gravity, m/s^2
rhog = 1.17                 # gas (air) density, kg/m^3

#------------------------------------------------------------------------------
# Expanded bed calculations

# bed expansion factor for dbed < 0.0635 m, Eq 14.7, pg 318 [Babu1978]
fbexp = 1 + (1.032*(u-umf)**0.57*rhog**0.083)/(rhop**0.166*umf**0.063*dbed**0.445)

# expanded bed height, zexp (m)
zexp = zmf*fbexp

# expanded bed void fraction, Eq 14.18, pg 320
eexp = 1 - (1-emf)/fbexp

#------------------------------------------------------------------------------
# Print results

print('')
print('zmf (m) =', zmf)
print('emf =', emf)
print('fbexp = ', fbexp)
print('zexp (m) =', zexp)
print('eexp =', eexp)