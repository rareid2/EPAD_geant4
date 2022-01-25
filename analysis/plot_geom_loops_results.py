import numpy as np
from plot_settings import *
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.patches import Ellipse, Circle 
import matplotlib.colors as colors
import seaborn as sns 

# ------------------------plotting function!-----------------------------
def plot_it(resolutions_list, labels, mt, color):

    sns.set_palette("Paired")
    
    # fix figure size
    fig, ax = plt.subplots(figsize=(6,3))


    for resolutions,label in zip(resolutions_list,labels):
        plt.plot(mask_distance,resolutions,label=label,linestyle='dashed', marker='s')
    
    plt.legend()
    plt.ylim([0,8])

    plt.xlabel('distance between mask and detector [cm]')
    plt.ylabel('angular resolution [deg]')
    plt.title(str(mt)+'um',color=color)

    abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"
    fname_save = 'results/fall21_results/ptsrcdiff/example_scatter.png'
    fname_save = 'results/fall21_results/ptsrcdiff/'+str(mt)+'_exp.png'

    fig.tight_layout(pad=0.3)
    plt.savefig(fname_save)
    plt.clf()


# parameter ranges
mask_thickness = np.linspace(400,3300,3) #um
mask_distance = np.linspace(0.5,6,5) #cm
mask_rank = [61,101,139,181]

# 8 means that it could not resolve 7!!

labels = ['61','101','139','181']
colors = ['#FFDC5E', '#DC493A', '#CB9CF2']


"""
# start with theoretical first
for mi,mt in enumerate(mask_thickness):
    th_resolutions = []
    labels = []
    for mr in mask_rank:
        pixel_size = round(63/mr,4)
        th_res = 2*np.rad2deg(np.arctan((pixel_size/2)/(np.array(mask_distance - (mt*1e-4))*10)))

        th_resolutions.append(th_res)
        labels.append(str(mr))


    plot_it(th_resolutions, labels, round(mt), colors[mi])
"""