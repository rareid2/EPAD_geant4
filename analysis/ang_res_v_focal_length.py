import numpy as np 
import matplotlib.pyplot as plt 
from plot_settings import * 
import os


# distance from CA to detector in mm
dist = np.linspace(0,90,100)
pixel_half = 0.055

mask_size = 28
detector_size = 14

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
fname = 'results/fall21_results/ang_res_fov_medipix_new.png'

fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()
"""

# new plan
mask_rank = np.linspace(61,167)
mask_element_size = 63/mask_rank # in mm

# stick w 15cm for now
geom_res = np.rad2deg(2*np.arctan((mask_element_size/2)/15))

plt.plot(mask_rank,geom_res)
fpath = os.getcwd()
fname = 'results/fall21_results/ang_res_v_fov/mask_size_v_res.png'
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()
"""