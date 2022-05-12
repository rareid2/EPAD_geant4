import numpy as np 
import matplotlib.pyplot as plt 
import os
from plot_settings import *
from fnc_find_theoretical_dist import findTheoreticalDist
from fnc_calc_angle_per_particle import calculateAnglePerParticle

# quick script to compare results to theoretical scattering MCS distribution (guassian)
thicknesses = np.logspace(0,3.5,10) # in um
colors = ["d9ed92","b5e48c","99d98c","76c893","52b69a","34a0a4","168aad","1a759f","1e6091","184e77"]
fig,ax = plt.subplots()


KEs = np.logspace(2,4,10)
for avg_KE,co in zip(KEs,colors):
    std_devs = []
    for det1_thickness_um in thicknesses:
        gap_in_cm = 3 #cm
        charge_nmbr = 1
        rest_mass_ME = 0.511 # kg
        # get data
        x_values, y_values, std_dev = findTheoreticalDist(det1_thickness_um, gap_in_cm, charge_nmbr, rest_mass_ME,avg_KE)

        #theta, theta_actual, avg_KE = calculateAnglePerParticle(gap_in_cm)
        #data_plt = theta

        std_devs.append(std_dev)

        # plot it!
        #ax.plot(x_values, y_values.pdf(x_values), '--',color='#'+co,label=str(round(det1_thickness_um)))

        #plt.hist(data_plt,density=True, bins=100,color=CB91_Blue,label='simulation results')
        #plt.yscale('log')
        #plt.ylim([0.0001,1])
    ax.plot(thicknesses, std_devs, '--',color='#'+co,label=str(round(avg_KE)))

plt.legend()
plt.margins(x=0)
plt.margins(y=0)
plt.ylim([0,180])
plt.xlim([1,100])

plt.xlabel('thickness of window [um]')
plt.ylabel('std of theoretical scattering dist. [deg]')
plt.title('electrons scattering from window')#  \n char length = ' + str(round(x/X0, 4)))
# save directory
folder_save = '/home/rileyannereid/workspace/geant4/EPAD_geant4/results/'
#fname = os.path.join(folder_save, 'updatedhighlanddistribution_results_'+str(int(round(avg_KE/1000)))+'MeV_'+str(det1_thickness_um)+'um.png')
fname = os.path.join(folder_save, 'updatedhighlanddistribution_results_Be.png')

plt.savefig(fname,dpi=800)
plt.close()

