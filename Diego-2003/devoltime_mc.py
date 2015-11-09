# Correlation for devolatilization time accounting for moisture content within
# a wood particle. Refer to Equations 2-7 in the paper for more details. 
# -- Reference (Diego2003) --
# De Diego, L. F., et al. "Effect of moisture content on devolatilization times
# of pine wood particles in a fluidized bed." Energy & fuels 17.2 (2003): 285-290.

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Parameters
# -----------------------------------------------------------------------------

d = 2       # diameter of particle, mm
n = 1.5     # exponent n
b = 1.7e-2  # factor for moisture content

# coefficient a at various temperatures
a1 = 1.69   # for 923 K
a2 = 1.38   # for 1023 K
a3 = 1.03   # for 1123 K
a4 = 0.89   # for 1223 K

# Calculate devolatilization time
#------------------------------------------------------------------------------

# range of moisture content from 0-50%
mc = np.linspace(0, 50) # moisture content, %

# general equation for devolatilization time tv = a*(d)^n per Eq. 2
# tv = tvo*(1+b*mc) per Eq. 5
tv1 = a1*d**n*(1 + b*mc)
tv2 = a2*d**n*(1 + b*mc)
tv3 = a3*d**n*(1 + b*mc)
tv4 = a4*d**n*(1 + b*mc)

# Plot results
#------------------------------------------------------------------------------

# this figure should have similar trends as Fig. 4 in Diego2003
py.figure(1)
py.plot(mc, tv1, lw=2, label='923 K')
py.plot(mc, tv2, lw=2, label='1023 K')
py.plot(mc, tv3, lw=2, label='1123 K')
py.plot(mc, tv4, lw=2, label='1223 K')
py.xlabel('Moisture Content (%)')
py.ylabel('tv (s)')
py.legend(loc='best', numpoints=1, fontsize=12)
py.title('d = {} mm particle size'.format(d))
py.grid()
py.show()

# plot the (a) coefficients vs temperature
T = [923, 1023, 1123, 1223]
a = [a1, a2 ,a3 ,a4]

py.figure(2)
py.scatter(T, a, s=60, c='red')
py.ylabel('Coefficient  a')
py.xlabel('Temperature (K)')
py.title('d = {} mm particle size'.format(d))
py.grid()
py.show()