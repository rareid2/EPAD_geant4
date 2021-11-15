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
colormap = plt.cm.YlGnBu
Ncolors=8
Ncolors = min(colormap.N,Ncolors)
mapcolors = [colormap(int(x*colormap.N/Ncolors)) for x in range(Ncolors)]
mapcolors = mapcolors[1:]
mapcolors.reverse()

# set energy from simulation 
all_energies = ['10','6','3p5','2','1','0p750','0p450']
ene_txt = ['10','6','3.5','2','1','0.75','0.45']
colors = ['palevioletred','mediumslateblue','lavender','indianred','lightskyblue','palegreen','cornflowerblue']
thicknesses = ['10','20','60','170']
txt_pos = [4.5,6.7,10.15,15.25,21.07,22.65,23.82]
rots = [3,6,11,17,11,8,2]
thicknesses_num = [10,20,60,170]

for enp,c,ent,ypos,ro in zip(all_energies,mapcolors,ene_txt,txt_pos,rots):
    save_std = []
    #rgba_color = cm.ylgnbl(norm(float(ene_txt)),bytes=True)

    for th in thicknesses:        
        fname = 'hits_'+enp+'_'+th+'.csv'
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
        theta_std = np.std(theta)

        save_std.append(theta_std)


    plt.plot(np.array(thicknesses_num), np.array(save_std),'--o',color=c)
    ax.text(125,ypos,ent+ ' MeV',fontsize='x-small',color=c,rotation=ro)


#plt.xlim([0, 200])
plt.ylabel('1-$\sigma$ [$^\circ$]')
plt.xlabel('detector thickness [$\mu$m]')
#plt.legend(loc='upper right')
plt.title('1-$\sigma$ uncertainty from mcs')
folder_save = '/home/rileyannereid/workspace/geant4/EPAD_geant4/results/fall21_results/reproducing_distribution'
fname = os.path.join(folder_save, 'one_sigma_v_thick.png')
plt.savefig(fname,dpi=300)