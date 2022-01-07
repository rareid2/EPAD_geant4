import numpy as np 
import os
import pandas as pd 
import scipy.stats
from scipy.signal import convolve2d as conv2
from scipy.fft import fft2,ifft
from skimage.transform import resize
from scipy.ndimage import zoom 

# function for reading hits
from fnc_get_det1_hits import getDet1Hits

# function for decoding matrix
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/rileyannereid/workspace/geant4/CA_designs')
from util_fncs import makeMURA, make_mosaic_MURA

# plotting
from plot_settings import *
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.patches import Ellipse, Circle 
import matplotlib.colors as colors

# ------------------------- ------------------------- ------------------------- -------------------------
# some helper functions

# Modified range to iterate over floats (i.e. 10.5 degrees, etc.)
def frange(start, stop, step):
     i = start
     while i < stop:
         yield i
         i += step

def shift(m, hs, vs):
    '''
    m: input image
    hs: horizontal shift
    vs: vertical shift
    '''
        
    # Get original image size
    rm,cm = np.shape(m);
    
    # Shift each quadrant by amount [hs, vs]
    m = np.block([[m[rm-vs:rm, cm-hs:cm], m[rm-vs:rm, 0:cm-hs]], 
                [m[0:rm-vs,cm-hs:cm], m[0:rm-vs,0:cm-hs]]])
        
    return m;

def fft_conv(rawIm, Dec):
    
    # scipy.ndimage.zoom used here
    resizedIm = zoom(rawIm, len(Dec)/len(rawIm));
        
    # Fourier space multiplication
    Image = np.real(np.fft.ifft2( np.fft.fft2(resizedIm) * np.fft.fft2(np.rot90(Dec, 2)) ));
    
    # Set minimum value to 0
    Image += np.abs(np.min(Image));
                
    # Shift to by half of image length after convolution
    return shift(Image, 16, 16);

# ------------------------- ------------------------- ------------------------- -------------------------
# settings 
abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"
fname_save = 'results/fall21_results/circular_beam/'
fname = 'data/hits.csv'

# get positions
xxes = []
yxes = []
fname_path = os.path.join(abs_path, fname)

# first get the x and z displacement
posX, posZ, energies = getDet1Hits(fname_path)

energy_limit_kev = 1

low_energy_electrons = 0
for x,z,ene in zip(posX,posZ,energies):
    if ene > energy_limit_kev:
        xxes.append(x)
        yxes.append(z)
    else:
        low_energy_electrons+=1

# plot the raw image
heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=63)
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
# remove outer 1mm on each side to match mask size for now
heatmap = heatmap[1:-1,1:-1]
plt.imshow(heatmap.T, extent=extent, origin='lower',cmap='turbo')
plt.colorbar(label='number of particles')
fname = fname_save + 'hits_raw.png'
fig_name = os.path.join(abs_path, fname)
plt.savefig(fig_name)
plt.close()

# ------------------------- ------------------------- ------------------------- -------------------------
# deconvolve the image
boxdim = 0.01 # cm -- just divide by 10 as well
nElements = 61 
holes_inv = True 
generate_files = False

mask, decode = makeMURA(nElements,boxdim,holes_inv,generate_files)

#mode='same'
#result_image = conv2(decode, heatmap.T, mode)
result_image = fft_conv(heatmap,decode)

c1 = plt.imshow(np.abs(result_image),cmap='turbo')
plt.colorbar(c1)
fname = fname_save + 'deconvolved.png'
fig_name = os.path.join(abs_path, fname)
plt.savefig(fig_name)
plt.close()
plt.clf()
# ------------------------- ------------------------- ------------------------- -------------------------