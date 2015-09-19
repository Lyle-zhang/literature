"""
Example of plotting data from a CSV file for Pyle 1984 data.
"""

# libraries and packages
import numpy as np
import matplotlib.pyplot as py

py.close('all')

# grab data from csv file
time, density = np.loadtxt('Fig6conv.csv', delimiter=',', unpack=True)

# plot the data
py.figure(1)
py.plot(time, density, 'og')
py.xlabel('Time (s)')
py.ylabel('Density (kg/m^3)')
py.grid()
py.show()