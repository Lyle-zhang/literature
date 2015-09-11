"""
Read in data from CSV file then plot.
"""

# libraries and packages
import numpy as np
import matplotlib.pyplot as py

py.close('all')

# parameters
Ti = 293    # initial temperature, K
Tf = 773    # ambient temperature, K

# grab data from csv file
time1, phi = np.loadtxt('Fig5b_phi.csv', delimiter=',', unpack=True)
temp = phi*(Ti-Tf)+Tf

time2, weight = np.loadtxt('Fig5b_weight.csv', delimiter=',', unpack=True)

# plot the data
py.figure(1)
py.plot(time1, phi, 'ob')
py.xlabel('Time (min)')
py.ylabel('Normalized Temperature (-)')
py.grid()

py.figure(2)
py.plot(time1, temp, 'or')
py.xlabel('Time (min)')
py.ylabel('Temperature (K)')
py.grid()

py.figure(3)
py.plot(time2, weight, 'og')
py.xlabel('Time (min)')
py.ylabel('Residual Weight Fraction (-)')
py.grid()

# -----------------------------------------------------------------------------

# parameters

Ti = 293
Tinf = 673

# grab data from csv files
t1, temp1 = np.loadtxt('Fig7center.csv', delimiter=',', unpack=True)
t2, temp2 = np.loadtxt('Fig7mid.csv', delimiter=',', unpack=True)
t3, temp3 = np.loadtxt('Fig7surf.csv', delimiter=',', unpack=True)

# plot the data
py.figure(4)
py.plot(t1, temp1+273, 'og')
py.plot(t2, temp2+273, 'ob')
py.plot(t3, temp3+273, 'or')
py.xlabel('Time (min)')
py.ylabel('Temperature (K)')
py.grid()

py.show()
