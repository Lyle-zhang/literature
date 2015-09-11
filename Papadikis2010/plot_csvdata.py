"""
Read in data from CSV file then plot.
"""

# libraries and packages
import numpy as np
import matplotlib.pyplot as py

py.close('all')

# grab data from csv file
time, density = np.loadtxt('Fig5_density350.csv', delimiter=',', unpack=True)

# plot the data
py.figure(1)
py.plot(time, density, 'ob')
py.xlabel('Time (s)')
py.ylabel('Density (kg/m^3)')
py.grid()
py.show()