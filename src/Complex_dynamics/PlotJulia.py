import numpy as np
import matplotlib.pyplot as plt

z = np.loadtxt("JuliaSet.dat")

plt.imshow(z,extent=[-2,2,-2,2])
plt.colorbar()
plt.show()
