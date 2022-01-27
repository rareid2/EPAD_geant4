import numpy as np 
import os
import pandas as pd 
import scipy.stats
from scipy.signal import convolve2d as conv2
from scipy.fft import fft2,ifft
from skimage.transform import resize
from scipy.ndimage import zoom 
import scipy.signal

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

    fig_name = os.path.join(abs_path, fname)
    plt.savefig(fig_name,bbox_inches='tight')
    plt.close()

def plot_peak(signal, fname):
    # convert the x axis???

    plt.plot(signal)

    # peak value
    half_val = (np.max(signal) - np.mean(signal[0:len(signal)//4]))//2
    plt.hlines(half_val + np.mean(signal[0:len(signal)//4]),xmin=0,xmax=len(signal),linestyles='--',colors='lightsalmon')

    fig_name = os.path.join(abs_path, fname)
    plt.savefig(fig_name,bbox_inches='tight')
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
    #resizedIm = rawIm
    #Dec = zoom(Dec, len(rawIm)/len(Dec));

    #fname = 'results/fall21_results/deconvolution_steps/resizedIm.png'
    #plot_step(resizedIm, vmax=150, fname=fname, label='# particles')

    Dec = np.fliplr(np.flipud(Dec))

    # Fourier space multiplication
    Image = np.real(np.fft.ifft2(np.fft.fft2(resizedIm) * np.fft.fft2(Dec) ));
    
    #fname = 'results/fall21_results/deconvolution_steps/postifft.png'
    #plot_step(Image, vmax=np.amax(Image), fname=fname, label='')

    # Set minimum value to 0
    Image += np.abs(np.min(Image));

    #fname = 'results/fall21_results/deconvolution_steps/minvalshift.png'
    #plot_step(Image, vmax=np.amax(Image), fname=fname, label='')
                
    # Shift to by half of image length after convolution
    return shift(Image, len(Dec)//2,len(Dec)//2)

# ------------------------- ------------------------- ------------------------- ------------------------- -------------------------------
# settings

def dec_plt(fname,uncertainty,nElements,boxdim,ff):

    abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"

    fname_save = 'results/fall21_results/ptsrcdiff/'+ ff
    #fname = 'data/hits.csv'

    #uncertainty = False
    #nElements = 11
    #boxdim = 3.0

    holes_inv = True
    generate_files = False

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
            if uncertainty:
                xxes.append(x + random.uniform(-1, 1)/10)
                yxes.append(z + random.uniform(-1, 1)/10)
            else:
                xxes.append(x)
                yxes.append(z)

        else:
            low_energy_electrons+=1

    # ------------------------------- ------------------------------- ------------------------------- -------------------------------
    heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=nElements)
    # remove outer 1mm on each side to match mask size for now
    #heatmap = heatmap[1:-1,1:-1]

    # plot the raw image
    #fname = fname_save + 'hits_raw.png'
    #plot_step(heatmap, vmax=np.amax(heatmap), fname=fname,label='# particles')

    # get decode array
    #mask, decode = makeMURA(nElements,boxdim,holes_inv,generate_files)
    mask, decode = make_mosaic_MURA(nElements,boxdim,holes=True,generate_files=False)
    decode = get_decoder_MURA(mask,nElements,holes_inv=False)

    # flip the heatmap over both axes bc point hole
    rawIm = np.fliplr(np.flipud(heatmap))

    # plot the flipped image
    #fname = fname_save + 'hits_flipped.png'
    #plot_step(rawIm, vmax=150, fname=fname,label='# particles')

    # reflect bc correlation needs to equal convolution
    rawIm = np.fliplr(rawIm)

    # plot the flipped image
    #fname = fname_save + 'hits_flipped_conv.png'
    #plot_step(rawIm, vmax=150, fname=fname,label='# particles')

    # add in poisson noise
    #noise_mask = np.random.poisson(rawIm)
    noise_mask = rawIm

    # deconvolve
    result_image = fft_conv(noise_mask,decode)
    snr = np.amax(result_image)/np.std(result_image)

    #print('SNR is: ', np.amax(result_image)/np.std(result_image))

    # plot final result
    #fname = fname_save + 'hits_deconvolve.png'
    #plot_step(result_image, vmax=np.amax(result_image), fname=fname,label='')

    # line plot of the diagonal -- must be flipped left and right but i have no clue why
    fname = fname_save + 'hits_deconvolve_peaks.png'
    plot_peak(np.diagonal(np.fliplr(result_image)),fname)

    return snr



"""
mask_thickness = np.linspace(400,3300,3) #um
mask_distance = np.linspace(0.5,6,5) #cm
mask_rank = [61,101,139,181]

mask_distance = ['0.5','1.875','3.25','4.625','6.0']

mr = str(181)
for mt in mask_thickness:
    resolutions = []
    mtt = str(round(mt))
    for md in mask_distance:
        fname = 'data/hits'+mr+'_'+mtt+'_'+md+'_0.csv'

        uncertainty = True
        nElements = int(mr)
        boxdim = round(63/nElements,4) #0.00649 # cm -- just divide by 10 as well
        boxdim = round(boxdim/100,6)
        holes_inv = True # keep false for mosaic
        generate_files = False

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
                if uncertainty:
                    xxes.append(x + random.uniform(-1, 1)/10)
                    yxes.append(z + random.uniform(-1, 1)/10)
                else:
                    xxes.append(x)# + random.uniform(-1, 1)/10)
                    yxes.append(z)# + random.uniform(-1, 1)/10)

            else:
                low_energy_electrons+=1

        # ------------------------------- ------------------------------- ------------------------------- -------------------------------
        heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=nElements)
        # remove outer 1mm on each side to match mask size for now
        #heatmap = heatmap[1:-1,1:-1]

        # plot the raw image
        fname = fname_save + 'hits_raw.png'
        plot_step(heatmap, vmax=np.amax(heatmap), fname=fname,label='# particles')

        # get decode array
        mask, decode = makeMURA(nElements,boxdim,holes_inv,generate_files)
        #mask, decode = make_mosaic_MURA(nElements,boxdim,holes=True,generate_files=False)
        #decode = get_decoder_MURA(mask,nElements,holes_inv)
        #print(decode)

        # flip the heatmap over both axes bc point hole
        rawIm = np.fliplr(np.flipud(heatmap))

        # plot the flipped image
        #fname = fname_save + 'hits_flipped.png'
        #plot_step(rawIm, vmax=150, fname=fname,label='# particles')

        # reflect bc correlation needs to equal convolution
        rawIm = np.fliplr(rawIm)

        # plot the flipped image
        #fname = fname_save + 'hits_flipped_conv.png'
        #plot_step(rawIm, vmax=150, fname=fname,label='# particles')

        # add in poisson noise
        #noise_mask = np.random.poisson(rawIm)
        noise_mask = rawIm

        # deconvolve
        result_image = fft_conv(noise_mask,decode)
        #print('SNR is: ', np.amax(result_image)/np.std(result_image))

        resolutions.append(np.amax(result_image)/np.std(result_image))

        # plot final result
        fname = fname_save + 'hits_deconvolve.png'
        plot_step(result_image, vmax=np.amax(result_image), fname=fname,label='')

        # line plot of the diagonal -- must be flipped left and right but i have no clue why
        fname = fname_save + 'hits_deconvolve_peaks.png'
        plot_peak(np.diagonal(np.fliplr(result_image)),fname)

    print(resolutions)

        # ------------------------- ------------------------- ------------------------- -------------------------
"""
# SNR first estimate:
"""
n_p = [100,193,373,720,1389,2683,5179,10000,19307,37276,71969,138950,268270,517947,1000000,1580000]
snr_n = [5.66,6.43,7.45,7.02,9.01,8.97,9.03,12.15,14.85,27.53,33.19,44.28,48.25,51.05,52.27,53.05]

fig, ax = plt.subplots()
plt.semilogx(n_p,snr_n,'.')
plt.xlabel('# particles')
plt.ylabel(' est. SNR')
fig.tight_layout(pad=0.5)
plt.savefig(fname_save + 'snr.png')
plt.close()

#100 - 5.66
#193 - 6.43
#373 - 7.45
#720 - 7.02
#1389 - 9.01
#2683 - 8.97
#5179 - 9.03
#10000 - 12.15
#19307 - 14.85
#37276 - 27.53
#71969 - 33.19
#138950 - 44.28
#268270 - 48.25
#517947 - 51.05
#1000000 - 52.27
#1580000 - 53.05
"""