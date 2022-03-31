
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm # To set the colormap in the plot.


data = np.loadtxt('wave.dat')

nt = np.shape(data)[0]
Nx = 100
Ny = 100

datamin = -.5 #np.min(data)
datamax = .5 #np.max(data)

def SliderPlot(imin, imax):
    # Import the function to allow the interactive slider:
    from matplotlib.widgets import Slider
    # Create a plot window with two objects in it (the image and the slider)
    fig, ax = plt.subplots()
    # Define the initial representation of the field:
    #frame = data[0,:].reshape(100,100)
    imag = plt.imshow(np.arcsinh(data[0,:]).reshape(100,100), vmin = np.arcsinh(datamin), vmax = np.arcsinh(datamax), cmap=cm.PuOr, origin='lower', extent=(0,14,0,10))
    # Make space in the window for the slider:
    plt.subplots_adjust(left=0.25, bottom=0.25)
    # Make the box to hold the time slider:
    axtime = plt.axes([0.25, 0.1, 0.65, 0.03])
    # Define the time slider object:
    global timeslider
    timeslider = Slider(axtime, "Time", imin, imax, valinit=0., valfmt='%1.2f')
    # Note that it has to be defined "globally" or else the animation would stop
    # working as soon as this function is no longer being evaluated.
    #
    # Now, define how the plot is updated if I adjust the time slider:
    def update(val):
        timeval = timeslider.val
        imag.set_data(np.arcsinh(data[int(timeval),:].reshape(100,100)))
        fig.canvas.draw_idle()
    # Prepare to update the plot anytime the slider is moved:
    timeslider.on_changed(update)
    # And show it:
    plt.show()
    


SliderPlot(0,nt-1)
