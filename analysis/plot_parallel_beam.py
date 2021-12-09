import numpy as np 
import os
import pandas as pd 
import scipy.stats
from scipy.signal import convolve2d as conv2

# function for reading hits
from fnc_get_det1_hits import getDet1Hits

# function for decoding
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/rileyannereid/workspace/geant4/CA_designs')
from util_fncs import makeMURA, make_mosaic_MURA

# plotting
from plot_settings import *
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.patches import Ellipse, Circle 
import matplotlib.colors as colors

# get positions
xxes = []
yxes = []

# energy limit
energy_limit_kev = 1

# define filename

fname = 'data/hits_5mm_250um_13mosaic.csv'
fdirectory = os.getcwd()
fname_path = os.path.join(fdirectory, fname)

# first get the x and z displacement from SCATTERING
posX, posZ, energies = getDet1Hits(fname_path)

low_energy_electrons = 0
for x,z,ene in zip(posX,posZ,energies):
    if ene > energy_limit_kev:
        xxes.append(x)
        yxes.append(z)
    else:
        low_energy_electrons+=1

#print('# of low energy electrons= ',low_energy_electrons)
heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=63)

extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

# deconvolve the image (?)
gridSizeX = 13
boxdim = 0.01  
mask, decode = make_mosaic_MURA(gridSizeX, boxdim)

# to plot the deconvolved image
mode='same'
corr = conv2(heatmap.T, decode, mode)
c1 = plt.imshow(conv2(heatmap.T, decode, mode),cmap='turbo',vmin=-4000,vmax=20000)
plt.colorbar(c1)
fname = 'results/fall21_results/deconvolved_5mm_250um_13mosaic.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()
plt.clf()

# plot the raw image
plt.imshow(heatmap.T, extent=extent, origin='lower',cmap='turbo',vmin=0,vmax=200)
plt.colorbar(label='number of particles')


fname = 'results/fall21_results/det_hits_13mosaic_5mm_250um.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()

#print(np.amax(corr))

# finished - plot the results of the window variation simulation - need to try with Be as well..... (tomorrow?)
#thick = np.array([50,150,250])
thick = np.array([20,60,100])

distance = np.array([1,3,5])
thickthick, dd = np.meshgrid(thick, distance, sparse=True)
correlation = np.array([[19391,16479,14341],[19404,4286,2382],[5410,1763,1431]])
#correlation = np.array([[29143,19827,14720],[15964,5736,4069],[10557,3017,1504]])

h = plt.contourf(thick,distance,correlation,vmin=0,vmax=32000)
plt.colorbar(h,label='correlation factor')
h.set_clim(0,32000)
plt.xlabel('thickness of Al window in um')
plt.ylabel('distance from window to detector in mm')
#plt.axis('scaled')
fname = 'results/fall21_results/window_sim_results_Al.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()
