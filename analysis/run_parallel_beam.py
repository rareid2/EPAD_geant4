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
fname = 'data/hits.csv'
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

# show the plot
plt.scatter(xxes, yxes,s=3)

fname = 'results/fall21_results/det1_hits_61array.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
plt.close()