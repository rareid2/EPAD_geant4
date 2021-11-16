import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.patches import Ellipse, Circle 
import scipy.stats
import os
from fnc_get_det_hits import getDetHits
from plot_settings import *


# set energy from simulation 
xxes = []
yxes = []

fname_path = os.path.join('/home/rileyannereid/workspace/geant4/Geant4_electron_detector/analysis/data', fname)
print(fname_path)
# first get the x and z displacement from SCATTERING
detector_hits, deltaX_rm, deltaZ_rm, energies = getDetHits(fname_path)
plt.scatter(deltaX_rm,deltaZ_rm)
plt.show()
plt.close()