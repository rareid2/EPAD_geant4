import numpy as np 
import matplotlib.pyplot as plt 
import os
from plot_settings import *

# quick script to compare results to theoretical scattering MCS distribution (guassian)
det1_thickness_um = 100 #um
gap_in_cm = 3 #cm
charge_nmbr = 2
rest_mass_kg = ##
x_values, y_values = findTheoreticalDist(det1_thickness_um, gap_in_cm, charge_nmbr, rest_mass_kg)

# plot it!
plt.plot(x_values, y_values.pdf(x_values), '--',color=CB91_Violet,label='theoretical')
plt.hist(data_plt,density=True, bins=100,color=CB91_Blue,label='simulation results')
plt.legend()
plt.margins(x=0)
plt.margins(y=0)
plt.xlabel('deg')
plt.title(str(int(avg_KE))+'MeV alpha scattering through detector 1  \n char length = ' + str(round(x/X0, 4)))
# save directory
folder_save = '/home/rileyannereid/workspace/geant4/Geant4_electron_detector/results/fall21_results/reproducing_distribution'
fname = os.path.join(folder_save, 'distribution_results_alpha_'+str(int(avg_KE))+'MeV_'+str(thickness)+'um.png')
plt.savefig(fname)
plt.close()