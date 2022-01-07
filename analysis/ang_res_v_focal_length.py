import numpy as np 
import matplotlib.pyplot as plt 
from plot_settings import * 
import os

# distance from CA to detector in mm
dist = np.linspace(0,80,100)
pixel_half = 1.0

mask_size = 70
detector_size = 63

# angular resolution in deg
theta = 2*np.arctan(pixel_half/dist)
theta = np.rad2deg(theta)

# FOV
fov = 2*np.arctan(((mask_size-detector_size)/2)/dist)
fov = np.rad2deg(fov)

ax1 = plt.subplot()

plt.plot(dist,theta,'b--')
ax1.set_ylim([0, 20])
plt.xlabel('mm between CA and detector')
plt.ylabel('theoretical limit of geom. res. [deg]',color='b')

ax2 = ax1.twinx()
ax2.plot(dist,fov,'r-')
plt.ylabel('fov [deg]', color='r')

fpath = os.getcwd()
fname = 'results/fall21_results/ang_res_microndetector_big_tp.png'

fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()