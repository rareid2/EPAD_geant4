import numpy as np 
import matplotlib.pyplot as plt 
import scipy.stats
import os
from plot_settings import *

# just a quick script to plot the 1 sigma uncertainty in theta

# open the results file
afile = open('/home/rileyannereid/workspace/geant4/EPAD_geant4/data/results.txt')

th = []
th_m = []
un = []
for li,line in enumerate(afile): 
    if li ==0:
        pass
    else:
        lines = line.split(',')
        th.append(float(lines[1]))
        th_m.append(float(lines[3]))
        un.append(abs(float(lines[4])))

# plotting
th_d = np.array(th_m) - np.array(th)
un = np.array(un)
plt.plot(np.linspace(0,30,60),th_d)
plt.ylim([-20, 20])
plt.plot(np.linspace(0,30,60),np.zeros(len(th_d)),'r--')
# Shade the area between y1 and y2
plt.fill_between(np.linspace(0,30,60), -1*un+th_d, un+th_d,
                 facecolor="purple", # The fill color
                 color='purple',       # The outline color
                 alpha=0.2)          # Transparency of the fill

plt.savefig('/home/rileyannereid/workspace/geant4/EPAD_geant4/results/fall21_results/reproducing_proposal_results/10MeV_QBBC_sens_det.png')
