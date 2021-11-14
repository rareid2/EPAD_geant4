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
import brewer2mpl

gap_in_cm = 3.0 # cm gap 

ax = plt.subplot(111)

# set energy from simulation 
all_energies = ['10','6','3p5','2','1','0p750','0p450']
ene_txt = ['10','6','3.5','2','1','0.75','0.45']
colors = ['palevioletred','mediumslateblue','lavender','indianred','lightskyblue','palegreen','cornflowerblue']
thicknesses = [0.04, 0.04001, 0.0402, 0.0405, 0.041, 0.038, 0.038]
text_pos = [(1,0.23),(1,0.15),(4,0.09),(6,0.05),(11,0.022),(14,0.015),(20,0.01)]

for enp,c,th,tp, ent in zip(all_energies,colors, thicknesses, text_pos,ene_txt):

    fname = 'hits_'+enp+'.csv'
    fname_path = os.path.join('/home/rileyannereid/workspace/geant4/EPAD_geant4/data', fname)
    print(fname_path)
    # first get the x and z displacement from SCATTERING
    detector_hits, deltaX_rm, deltaZ_rm, energies = getDetHits(fname_path)
    
    x = deltaZ_rm
    
    # add in uncertainty from the position

    newx = []
    for dx in x:
        dp = random.uniform(-np.sqrt(2), np.sqrt(2))
        dx = (dp/10)+dx
        newx.append(dx)

    x = np.array(newx)

    # Find angles in degrees
    theta = np.rad2deg(np.arctan2(x, gap_in_cm))

    plt.hist(theta, bins=100, density=True, alpha=0.5,
         histtype='stepfilled', color=c,
         edgecolor='none',label=str(round(np.std(theta),1)) +'$^\circ$')

    # quick script to compare results to theoretical scattering MCS distribution (guassian)
    det1_thickness_um = th*10**3 #um
    gap_in_cm = 3 #cm
    charge_nmbr = 1
    rest_mass_ME = 0.511 # kg
    # get data
    x_values, y_values = findTheoreticalDist(det1_thickness_um, gap_in_cm, charge_nmbr, rest_mass_ME,fname_path)

    plt.plot(x_values, y_values.pdf(x_values), '--',color=c)
    plt.text(tp[0],tp[1],ent,color=c,fontsize=8,rotation=30)

plt.xlim([-70, 70])
plt.xlabel('theta [deg]')
plt.legend(loc='upper right')
plt.title('angular distribution of scattering through detector 1')
folder_save = '/home/rileyannereid/workspace/geant4/EPAD_geant4/results/fall21_results/reproducing_distribution'
fname = os.path.join(folder_save, 'all_scattering.png')
plt.savefig(fname,dpi=300)