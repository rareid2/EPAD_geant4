import numpy as np 
import matplotlib.pyplot as plt 
import os
from plot_settings import *
from fnc_find_theoretical_dist import findTheoreticalDist
from fnc_calc_angle_per_particle import calculateAnglePerParticle

# quick script to compare results to theoretical scattering MCS distribution (guassian)
det1_thickness_um = 40 #um
gap_in_cm = 3 #cm
charge_nmbr = 1
rest_mass_ME = 0.511 # kg
# get data
x_values, y_values = findTheoreticalDist(det1_thickness_um, gap_in_cm, charge_nmbr, rest_mass_ME)
theta, theta_actual, avg_KE = calculateAnglePerParticle(gap_in_cm)
data_plt = theta

# plot it!
plt.plot(x_values, y_values.pdf(x_values), '--',color=CB91_Violet,label='theoretical')
plt.hist(data_plt,density=True, bins=100,color=CB91_Blue,label='simulation results')
#plt.yscale('log')
#plt.ylim([0.0001,1])
plt.legend()
plt.margins(x=0)
plt.margins(y=0)
plt.xlabel('deg')
plt.title(str(int(round(avg_KE/1000)))+'MeV electrons scattering through detector 1')#  \n char length = ' + str(round(x/X0, 4)))
# save directory
folder_save = '/home/rileyannereid/workspace/geant4/EPAD_geant4/results/fall21_results/reproducing_distribution'
fname = os.path.join(folder_save, 'updatedhighlanddistribution_results_'+str(int(round(avg_KE/1000)))+'MeV_'+str(det1_thickness_um)+'um.png')
plt.savefig(fname)
plt.close()