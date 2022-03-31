
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("wave.dat")

#print(np.shape(data))

nt, nx = np.shape(data)

#for i in range(nt):
#    plt.plot(data[i,:])
#    plt.ylim(-.2,1.2)
#    plt.show()


plt.ion()
curve, = plt.plot(data[0,:])
plt.ylim(-1.3, 1.3)
#plt.show()

#plt.ion()

for i in range(nt):
    curve.set_ydata(data[i,:])
    plt.draw()
    plt.pause(0.001)

plt.ioff()
