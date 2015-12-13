import numpy as np

# load the particle position data
# resulting data:
# pos[ <number timesteps> , <number particle> , 2 (spatial dimensions) ]
print 'Loading particle data...'
execfile('VTK/timestep_all_particles.py')

# print some information
print 'Done.'
print 'Number of timesteps: \t' + str(pos.shape[0])
print 'Number of particles: \t' + str(pos.shape[1])

# plot
import matplotlib.pyplot as plt
import time
from math import floor, ceil

xmin = floor(np.nanmin(pos[:,:,0]))
xmax = ceil(np.nanmax(pos[:,:,0]))
ymin = floor(np.nanmin(pos[:,:,1]))
ymax = ceil(np.nanmax(pos[:,:,1]))

for ii in range(pos.shape[0]):
    plt.plot(pos[ii,:,0], pos[ii,:,1], 'bo', hold=False)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.draw()
    plt.show(block=False)
    time.sleep(0.02)

plt.show(block=True)