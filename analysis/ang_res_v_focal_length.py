import numpy as np 
import matplotlib.pyplot as plt 
from plot_settings import * 
import os

# distance from CA to detector in mm
dist = np.linspace(0,100,100)
pixel_half = 0.5

# angular resolution in deg
theta = 2*np.arctan(pixel_half/dist)
theta = np.rad2deg(theta)

plt.plot(dist,theta,'--')
plt.xlabel('mm between CA and detector')
plt.ylabel('theoretical limit of detector resolution [deg]')

fpath = os.getcwd()
fname = 'results/fall21_results/ang_res.png'
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()

# angular resolution in deg
dist = np.linspace(10,50,100)
theta = 2*np.arctan(pixel_half/dist)
theta = np.rad2deg(theta)

plt.plot(dist,theta,'--')
plt.xlabel('mm between CA and detector')
plt.ylabel('theoretical limit of detector resolution [deg]')

fname = 'results/fall21_results/ang_res_zoom.png'
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()