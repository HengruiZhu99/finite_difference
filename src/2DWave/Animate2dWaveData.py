
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm


data = np.loadtxt('wave.dat')

nt = np.shape(data)[0]

Nx = 100
Ny = 100

plt.ion()
graph = plt.imshow(data[0,:].reshape(Nx,Ny), origin='lower', cmap = cm.PuOr, extent=(0,14,0,10), vmin=-.5, vmax=.5)
plt.draw()

for i in range(nt):
    graph.set_data(data[i,:].reshape(Nx, Ny))
    plt.draw()
    plt.pause(0.01)

plt.ioff()
