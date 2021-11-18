import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.patches import Ellipse, Circle 
import scipy.stats
import os
from fnc_get_det1_hits import getDet1Hits
from plot_settings import *

# get positions
xxes = []
yxes = []

# energy limit
energy_limit_kev = 400

# define filename
fname = 'data/hits_500keV_100e3_windowbehind.csv'
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

plt.clf()
plt.imshow(heatmap.T, extent=extent, origin='lower',cmap='gist_gray')
fname = 'results/fall21_results/det1_hits_61array.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()