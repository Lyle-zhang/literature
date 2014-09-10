# Correlation for devolatilization time as tv
# where tv = time for 95% conversion
# Blasi2003, Eq. 2, based on Figure 12

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Parameters
# -----------------------------------------------------------------------------

Tr1 = 700   # reactor temperature, K
Tr2 = 800
Tr3 = 900
Tr4 = 1000
d = np.linspace(0.5, 11)  # particle diameter, mm

# Calculate devolatilization time, tv
# -----------------------------------------------------------------------------

# tv (devolatilization time) from Eq. 2 in Blasi2003
# it is the time when conversion is about 95%
tv1 = 0.8*np.exp(1525/Tr1)*(d**1.2) # devol. time for Tr1
tv2 = 0.8*np.exp(1525/Tr2)*(d**1.2) # devol. time for Tr2
tv3 = 0.8*np.exp(1525/Tr3)*(d**1.2) # devol. time for Tr2
tv4 = 0.8*np.exp(1525/Tr4)*(d**1.2) # devol. time for Tr2

# Plot results
# -----------------------------------------------------------------------------

py.figure(1)
py.plot(d, tv1)
py.text(max(d), max(tv1)+1, '700K')
py.plot(d, tv2)
py.text(max(d), max(tv2)+1, '800K')
py.plot(d, tv3)
py.text(max(d), max(tv3)+1, '900K')
py.plot(d, tv4)
py.text(max(d), max(tv4)+1, '1000K')
py.xlabel('d (mm)')
py.ylabel('tv (s)', rotation='horizontal', horizontalalignment='right')
py.title('95% conversion time')
py.grid()
py.show()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Parameters
# -----------------------------------------------------------------------------

d1 = 0.5  # particle diameter, mm
d2 = 1.0
d3 = 2
d4 = 4
Tr = np.linspace(625, 1125) # reactor temperature, K

# Calculate devolatilization time, tv
# -----------------------------------------------------------------------------

tv1 = 0.8*np.exp(1525/Tr)*(d1**1.2)
tv2 = 0.8*np.exp(1525/Tr)*(d2**1.2)
tv3 = 0.8*np.exp(1525/Tr)*(d3**1.2)
tv4 = 0.8*np.exp(1525/Tr)*(d4**1.2)

# Plot results
# -----------------------------------------------------------------------------

py.figure(2)
py.plot(Tr, tv1)
py.text(max(Tr)+5, min(tv1)-0.5, str(d1)+' mm')
py.plot(Tr, tv2)
py.text(max(Tr)+5, min(tv2)-0.5, str(d2)+' mm')
py.plot(Tr, tv3)
py.text(max(Tr)+5, min(tv3)-0.5, str(d3)+' mm')
py.plot(Tr, tv4)
py.text(max(Tr)+5, min(tv4)-0.5, str(d4)+' mm')
py.xlabel('Tr (K)')
py.ylabel('tv (s)', rotation='horizontal', horizontalalignment='right')
py.title('95% conversion time')
py.grid()
py.show()