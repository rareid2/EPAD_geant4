import numpy as np 
import os
import pandas as pd 
import scipy.stats
from scipy.signal import convolve2d as conv2
from scipy.fft import fft2,ifft
from skimage.transform import resize
from scipy.ndimage import zoom 
import scipy.signal
import matplotlib.pyplot as plt
import random 

# function for reading hits
from fnc_get_det1_hits import getDet1Hits

# function for decoding matrix
import sys
# insert at 1, 0 is the script path (or '' in REPL)
sys.path.insert(1, '/home/rileyannereid/workspace/geant4/CA_designs')
from util_fncs import makeMURA, make_mosaic_MURA, get_decoder_MURA

# plotting
from plot_settings import *
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.patches import Ellipse, Circle 
import matplotlib.colors as colors

# ------------------------- ------------------------- ------------------------- ------------------------- -------------------------------
# some helper functions

def plot_step(signal, vmax, fname, label):
    plt.imshow(signal.T, origin='lower',cmap='turbo',vmax=vmax)
    plt.colorbar(label=label)
    #plt.xlim([40,55])
    #plt.ylim([40,55])

    plt.savefig(fname,bbox_inches='tight')
    plt.close()

def plot_peak(signal, fname):
    plt.plot(signal)

    # peak value
    half_val = (np.max(signal) - np.mean(signal[0:len(signal)//4]))//2
    plt.hlines(half_val + np.mean(signal[0:len(signal)//4]),xmin=0,xmax=len(signal),linestyles='--',colors='lightsalmon')

    fname = fname + '.png'

    plt.savefig(fname,bbox_inches='tight')
    plt.close()

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
    hs+=1
    vs+=1
        
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
    Image = np.real(np.fft.ifft2(np.fft.fft2(resizedIm) * np.fft.fft2(Dec) ));
    
    # Set minimum value to 0
    Image += np.abs(np.min(Image));
                
    # Shift to by half of image length after convolution
    return shift(Image, len(Dec)//2,len(Dec)//2)

# ------------------------- ------------------------- ------------------------- ------------------------- -------------------------------
# settings

# fname = where is the data
# uncertainty = add uncertainty? 

def dec_plt(fname,uncertainty,nElements,boxdim,ff):
    
    abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"

    fname_save = 'results/fall21_results/ptsrcdiff/'+ ff

    # get positions
    xxes = []
    yxes = []
    fname_path = os.path.join(abs_path, fname)

    # first get the x and y displacement in cm
    posX, posY, energies = getDet1Hits(fname_path)

    energy_limit_kev = 1
    low_energy_electrons = 0

    # filter so that each pixel is a bin (remove outer edges)
    bins = 63 / boxdim
    bins_extra = bins - np.floor(bins)
    
    # remove the outer 
    extra = 2 # convert from mm to cm
    shift = 8.4/2
    x_out = []
    y_out = []

    for x,y,ene in zip(posX,posY,energies):
        # shift origin to lower left
        x += shift
        y += shift
        if ene > energy_limit_kev:
            if x < extra or y < extra or x > (8.4 - extra) or y > (8.4 - extra): 
                #print('oops')
                x_out.append(x)
                y_out.append(y)
            else:
                if uncertainty:
                    xxes.append(x + random.uniform(-1, 1)/10)
                    yxes.append(y + random.uniform(-1, 1)/10)
                else:
                    xxes.append(x)
                    yxes.append(y)

        else:
            low_energy_electrons+=1

    # ------------------------------- ------------------------------- ------------------------------- -------------------------------
    heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=11*4)

    fname_step = fname_save + '_raw.png'
    plot_step(heatmap,np.amax(heatmap),fname_step,label='# particles')

    # first get the mask to use in 
    mask, decode = make_mosaic_MURA(nElements,boxdim,holes=False,generate_files=False)

    fname_step = fname_save + '_mask1.png'
    
    plot_step(mask,np.amax(mask),fname_step,label='# particles')

    decode = get_decoder_MURA(mask,nElements,holes_inv=False)

    fname_step = fname_save + '_mask.png'
    
    plot_step(decode,np.amax(decode),fname_step,label='# particles')

    # flip the heatmap over both axes bc point hole
    rawIm = np.fliplr(np.flipud(heatmap))

    # reflect bc correlation needs to equal convolution
    rawIm = np.fliplr(rawIm)

    # deconvolve
    result_image = fft_conv(rawIm,decode)
    fname_step = fname_save + '_dc.png'
    plot_step(result_image,np.amax(result_image),fname_step,label='# particles')

    snr = np.amax(result_image)/np.std(result_image)

    # line plot of the diagonal -- must be flipped left and right but i have no clue why
    plot_peak(np.diagonal(np.fliplr(result_image)),fname_save)

    return snr

#dec_plt('data/hits.csv',False, 11, 3.0, 'test11')