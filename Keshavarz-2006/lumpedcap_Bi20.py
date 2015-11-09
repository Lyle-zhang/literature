# Improved lumped capacitance method for Bi <= 20
# from article Keshavarz2006

import numpy as np
import matplotlib.pyplot as py

#py.close('all')

# parameters

h = 350     # convection heat transfer coefficient, W/m^2*K
k = 0.12    # thermal conductivity, W/m*K
rho = 540   # density, kg/m^3
c = 2000    # heat capacity, J/kg*K
d = 700e-6  # particle diameter, m (e-6 for microns)

Ti = 300    # initial temperature, K
Tinf = 773  # ambient temperature, K

# calculations

t = np.linspace(0, 3, 100)   # time range, s

A = np.pi*(d**2)        # surface area sphere pi*D^2, m^2
V = np.pi*(d**3)/6      # volume of sphere pi*D^3/6, m^3

Lc = V/A                # characteristic length, m
Bi = (h*Lc)/k           # Biot number, ~

alpha = k/(rho*c)       # thermal diffusivity, m^2/s
Fo = (alpha*t)/(Lc**2)  # Fourier number, ~

m = 2                   # 0=slab, 1=cylinder, 2=sphere

tm = ((m+1)/(m+3))*Bi+1
phi = np.exp(-(1/tm)*Bi*Fo)     # dimensionless temperature, Eq 22
T = Tinf+(Ti-Tinf)*phi

# plot

py.figure(1)
py.plot(t, T, linewidth=2, label='lump')
py.axhline(y=773, color='k', linestyle='--')
py.ylim(200, 800)
py.xlim(0, 3)
py.title('Lumped Capacitance Method for Sphere\n for Bi <= 20, Keshavarz2006')
py.ylabel('Volume Avg. Temperature (K)')
py.xlabel('Time (s)')
py.legend(loc='best', numpoints=1)
py.grid()
py.show()