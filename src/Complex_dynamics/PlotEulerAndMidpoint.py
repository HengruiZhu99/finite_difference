
import numpy as np
import matplotlib.pyplot as plt


dataE = np.loadtxt("EulerStepping.dat")
dataM = np.loadtxt("MidpointStepping.dat")

tE = dataE[:,0]
uE = dataE[:,1]

tM = dataM[:,0]
uM = dataM[:,1]

plt.plot(tE, uE, label="Euler's method")
plt.plot(tM, uM, label="Midpoint method")
plt.show()
