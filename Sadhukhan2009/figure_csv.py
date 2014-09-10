# use Python 3 print function and division
from __future__ import print_function
from __future__ import division

# libraries and packages
import numpy as np
import matplotlib.pyplot as py

# grab data from csv file
m = np.loadtxt('Fig1_Mcylinder.csv', delimiter=',')

# plot the data
py.figure(1)
py.plot(m[:,0], m[:,1], 'ob')
py.xlabel('Time (s)')
py.ylabel('Mass Fraction (-)')
py.grid()
py.show()