import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.patches import Ellipse, Circle 
import scipy.stats
import os
from fnc_get_det_hits import getDetHits
from plot_settings import *

import random
import decimal
import brewer2mpl


def eigsorted(cov):
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:,order]

#truly unsure whats up here
def add_contours(x,y,color_contour,ax,enp):
    cov = np.cov(x, y)
    vals, vecs = eigsorted(cov)
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    
    for j in range(1,2):
        w, h = 2 * nstd * np.sqrt(vals)*j
        ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                    width=w, height=h,
                    angle=theta, color=color_contour)
        ell.set_facecolor('none')
        ax.add_artist(ell)
    print(enp)
    plt.text(0,(w/2)+0.02,enp,fontsize=8,color=color_contour)


# set energy from simulation 
xxes = []
yxes = []

for enp in np.logspace(2,4,10):
    energy = int(round(enp))

    enp = int(round(enp))

    fname = 'hits_'+str(enp)+'.csv'
    fname_path = os.path.join('/home/rileyannereid/workspace/geant4/Geant4_electron_detector/analysis/data', fname)
    print(fname_path)
    # first get the x and z displacement from SCATTERING
    detector_hits, deltaX_rm, deltaZ_rm, energies = getDetHits(fname_path)

    # get the hits as well -- since we know this is for straight on this is easier
    # make a scatter plot of detector 1
    """
    ax = plt.subplot(111)
    plt.scatter(0,0)

    plt.xlim([-3.15,3.15])
    plt.ylim([-3.15,3.15])
    plt.xlabel('cm')
    plt.ylabel('cm')
    ax.set_aspect('equal')
    plt.title(str(energy)+' detector 1 hits')
    folder_save = '/home/rileyannereid/workspace/geant4/Geant4_electron_detector/results/fall21_results/reproducing_distribution'
    fname = os.path.join(folder_save, str(energy)+'_detector1_hits.png')
    plt.savefig(fname)
    plt.close()
    """

    ## ---- detector 2 ---------------
    ax = plt.subplot(111)

    # draw the standard deviation on here -- elliptical 
    x = deltaX_rm
    y = deltaZ_rm

    # add in uncertainty from the position
    newx = []
    newy = []
    for dx, dy in zip(x,y):
        dp = random.uniform(-np.sqrt(2), np.sqrt(2))
        dn = random.uniform(-np.sqrt(2), np.sqrt(2))

        dx = (dp/10)+dx
        dy = (dn/10)+dy

        newx.append(dx)
        newy.append(dy)

    x = np.array(newx)
    y = np.array(newy)

    nstd = 2
    
    xxes.append(x)
    yxes.append(y)



#plt.scatter(x, y,s=np.ones(len(x)))

# plot ratio
# Get "Set2" colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)

ax = plt.subplot(111)
energies = np.logspace(2,4,10)
colors = ['palevioletred','indianred','mediumslateblue', 'palegreen','lightskyblue','lavender','grey','cornflowerblue','olive','salmon']
for x,y,cc,en in zip(xxes,yxes,colors,energies):
    enp = str(int(round(en)))
    add_contours(x,y,cc,ax,enp)


plt.xlim([-3.15,3.15])
plt.ylim([-3.15,3.15])
plt.xlabel('cm')
plt.ylabel('cm')
ax.set_aspect('equal')

# add in some info
#plt.text(-2.5,2.5,'1-sigma uncertainty in theta: '+ str(np.round(np.rad2deg(np.arctan(np.std(x)/3)), 2)) + ' deg',fontsize=12,color='red')
#plt.text(-2.5,2,'1-sigma uncertainty in phi: '+ str(np.round(np.rad2deg(np.arctan(np.std(y)/3)), 2)) + ' deg',fontsize=12,color='red')

# now add in uncertainty from the position resolution ??? 
#circle1 = Circle((0, 0), 0.1, color='r')
#ax.add_patch(circle1)
#circle1.set_facecolor('none')

plt.title('detector 2 hits')
folder_save = '/home/rileyannereid/workspace/geant4/Geant4_electron_detector/results/fall21_results/reproducing_distribution'
fname = os.path.join(folder_save, 'all_detector2_hits.png')
plt.savefig(fname,dpi=300)

# angular width of the uncertainty in position resolution
print(str(np.round(np.rad2deg(np.arctan(0.1/3)), 2)))