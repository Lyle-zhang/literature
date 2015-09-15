"""
Example of plotting data from a CSV file for Sadhukhan 2009.
"""

# libraries and packages
import numpy as np
import matplotlib.pyplot as py

py.close('all')

# grab data from csv file
t1, Tcyl = np.loadtxt('Fig1_Tcylinder.csv', delimiter=',', unpack=True)
t2, Mcyl = np.loadtxt('Fig1_Mcylinder.csv', delimiter=',', unpack=True)


# plot the data
py.figure(1)
py.plot(t1, Tcyl+273, 'og', label='axis')
py.legend(loc='best', numpoints=1)
py.xlabel('Time (s)')
py.ylabel('Temperature (K)')

py.figure(2)
py.plot(t2, Mcyl, 'og', label='weight')
py.ylim([0, 1.1])
py.legend(loc='best', numpoints=1)
py.xlabel('Time (s)')
py.ylabel('Residual Weight Fraction (-)')

py.show()