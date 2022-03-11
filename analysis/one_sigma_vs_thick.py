import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.patches import Ellipse, Circle 
import scipy.stats
import os
from fnc_get_det_hits import getDetHits
from plot_settings import *
from fnc_find_theoretical_dist import findTheoreticalDist
import random
import decimal
import seaborn as sns 

def get_data(fname):

    file1 = open(fname, 'r')
    
    # Using readlines()
    Lines = file1.readlines()
    
    count = 0
    # Strips the newline character
    all_data = []
    for line in Lines:
        count += 1
        data = line.strip()
        data = data.split(',')
        clean_data = [float(d) for d in data]
        all_data.append(clean_data)

    file1.close()
    return all_data

fname = 'results/parameter_sweeps/two_det_3.0.txt'
my_data = get_data(fname)
fname = 'results/parameter_sweeps/two_det_0.5.txt'
my_data2 = get_data(fname)

all_data = []
for di,dd in zip(my_data,my_data2):
    data = []
    data.append(di)
    data.append(dd)
    all_data.append(data)


thicknesses = np.linspace(20,200,10)
energies = np.logspace(2,4,10)
energies = reversed(energies[5:])

coolors = ['#012a4a','#013a63','#01497c','#014f86','#2a6f97','#2c7da0','#468faf','#61a5c2','#89c2d9','#a9d6e5']
coolors = reversed(coolors)
fig, ax = plt.subplots()
sns.set_palette("Paired")
for resolutions,color in zip(all_data,coolors):
    ax.plot(thicknesses,resolutions[0],linestyle='--', marker='.',zorder=1,color=color)
    ax.plot(thicknesses,resolutions[1],linestyle='--', marker='.',zorder=1,color=color)
    ax.fill_between(thicknesses, resolutions[0], resolutions[1], alpha=0.5,color=color)

ene_text = '#700353'

ax.set_ylabel('uncertainty in incident polar angle [deg]')
ax.set_xlabel('thickness of front detector [um]')
ax.grid(which='major', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()
ax.set_axisbelow(True)


custom_loc = [4.3,7,11,16,20.75]
custom_rot = [5,10,13,15,10]
ii = 0
for ene in energies:
    ene_st = str(int(ene))
    plt.text(180, custom_loc[ii], ene_st + ' keV',rotation=custom_rot[ii],fontsize=10,color=ene_text)
    ii+=1
plt.text(170, 23.3, '100-775 keV',rotation=0,fontsize=10,color=ene_text)

plt.savefig('results/parameter_sweeps/two_det_results.png',dpi=800)
plt.close()