# Correlation for devolatilization time as tv where tv = time for 95% conversion
# Blasi2003, Eq. 2, Figure 12 example

import numpy as np
import matplotlib.pyplot as py

py.close('all')

# Parameters
# -----------------------------------------------------------------------------

Tr1 = 807
Tr2 = 1107
d = np.linspace(1, 11)

# Calculate devolatilization time, tv
# -----------------------------------------------------------------------------

tv1 = 0.8*np.exp(1525/Tr1)*(d**1.2)
tv2 = 0.8*np.exp(1525/Tr2)*(d**1.2)

# Plot results
# -----------------------------------------------------------------------------

py.figure(1)
py.plot(d, tv1)
py.text(max(d)-0.5, max(tv1)+1, '807K')
py.plot(d, tv2)
py.text(max(d)-0.5, max(tv2)+1, '1107K')
py.xlabel('d (mm)')
py.ylabel('tv (s)')
py.title('Blasi2003 - Figure 12')
py.grid()
py.show()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# Parameters
# -----------------------------------------------------------------------------

d1 = 4
d2 = 8
Tr = np.linspace(625, 1125)

# Calculate devolatilization time, tv
# -----------------------------------------------------------------------------

tv1 = 0.8*np.exp(1525/Tr)*(d1**1.2)
tv2 = 0.8*np.exp(1525/Tr)*(d2**1.2)

# Plot results
# -----------------------------------------------------------------------------

py.figure(2)
py.plot(Tr, tv1)
py.text(max(Tr)+0.5, min(tv1)+1, '4 mm')
py.plot(Tr, tv2)
py.text(max(Tr)+0.5, min(tv2)+1, '8 mm')
py.xlabel('Tr (K)')
py.ylabel('tv (s)')
py.title('Blasi2003 - Figure 12')
py.grid()
py.show()