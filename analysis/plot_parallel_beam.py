import numpy as np 
import os
import pandas as pd 
import scipy.stats
from scipy.signal import convolve2d as conv2

# function for reading hits
from fnc_get_det1_hits import getDet1Hits

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

fname = 'data/hits_win_front_retry.csv'
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

print('# of low energy electrons= ',low_energy_electrons)
heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=63)

extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

# deconvolve the image (?) 
decode_file = '/home/rileyannereid/workspace/geant4/CA_designs/decoding_matrix.txt'
decoder = pd.read_csv(decode_file, delimiter=',')

# reconstruct the array as 2D
decode_array = np.zeros((61,61))
count = 0
for i in range(61):
    for j in range(61):
        decode_array[i,j] = decoder.d[count]
        count+=1

# to plot the decoder
#c1=plt.scatter(decoder.j, decoder.i, c=decoder.d, marker='s')
#plt.colorbar(c1);

# to plot the deconvolved image
mode='same'
c1 = plt.imshow(conv2(heatmap.T, decode_array, mode),cmap='turbo',vmin=-6000,vmax=35000)
plt.colorbar(c1)

fname = 'results/fall21_results/deconvolved_window_infront_retry.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()
plt.clf()

# plot the raw image
plt.imshow(heatmap.T, extent=extent, origin='lower',cmap='turbo',vmin=0,vmax=50)
plt.colorbar(label='number of particles')


fname = 'results/fall21_results/det_hits_61array_windowfront_retry.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()
