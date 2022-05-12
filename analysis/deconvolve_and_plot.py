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
import gc 

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
    plt.clf()

def plot_peak(signal, fname, condition):
    # create an x axis
    x_ax = np.arange(0,len(signal))

    # plot the signal
    plt.plot(signal)

    # plot local mins and maxes
    b = (np.diff(np.sign(np.diff(signal))) > 0).nonzero()[0] + 1 
    c = (np.diff(np.sign(np.diff(signal))) < 0).nonzero()[0] + 1

    plt.scatter(x_ax[b],signal[b],color='b')
    plt.scatter(x_ax[c],signal[c],color='r')

    # define conditions for half and quarter separated
    half_val = (np.max(signal) - np.mean(signal[0:len(signal)//4]))//2
    half_val = half_val + np.mean(signal[0:len(signal)//4])
    quarter_val = (np.max(signal) - np.mean(signal[0:len(signal)//4]))//4
    quarter_val = 3*quarter_val + np.mean(signal[0:len(signal)//4])

    # plot them
    plt.hlines(half_val,xmin=0,xmax=len(signal),linestyles='--',colors='lightsalmon')
    plt.hlines(quarter_val,xmin=0,xmax=len(signal),linestyles='--',colors='lightsalmon')

    # find the height of the two peaks (if there are two peaks)
    local_maxes = signal[c]
    local_maxes.sort()  
    largest_peak = local_maxes[-1]
    second_largest_peak = local_maxes[-2]
    largest_local_min = np.max(signal[b])

    if condition=='half_val':
        # first, are there two peaks? 
        # peak here is defined as larger than half value
        if largest_peak > half_val and second_largest_peak > half_val:
            # is the largest local min below the condition
            if largest_local_min < half_val:
                resolution = True
            else:
                resolution = False
                #print('peaks are separated but not enough')
        else:
            resolution = False
            #print('peaks are not separated')

    elif condition=='quarter_val':
        if largest_peak > quarter_val and second_largest_peak > quarter_val:
            # is the largest local min below the condition
            if largest_local_min < quarter_val:
                resolution = True
            else:
                resolution = False
        else:
            resolution = False

    fname = fname + '.png'

    plt.savefig(fname,bbox_inches='tight')
    plt.close()
    plt.clf()

    return resolution

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

# TO DO 
# change how uncertainty is added in

def dec_plt(fname,uncertainty,nElements,boxdim,ff,ms):
    
    abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"

    fname_save = 'results/parameter_sweeps/'+ ff

    # needs to be set for decoding
    if nElements == 67 or nElements == 31:
        check =0
    else:
        check=1

    # get positions
    xxes = []
    yxes = []
    fname_path = os.path.join(abs_path, fname)

    # first get the x and y displacement in cm
    posX, posY, energies = getDet1Hits(fname_path)

    energy_limit_kev = 1
    low_energy_electrons = 0
    
    # remove the outer
    if nElements == 61:
        multiplier = 4
    else:
        multiplier = 8

    multiplier = 12
    #detector_sz = 1.364
    detector_sz = 4.422
    #multiplier = 30
    
    out_size = round((detector_sz-(nElements*boxdim/10))/2,3) # convert from mm to cm

    shift = detector_sz/2
    rr = []

    for x,y,ene in zip(posX,posY,energies):
        # shift origin to lower left
        x += shift
        y += shift
        if ene > energy_limit_kev:
            # confirm not out of bounds
            if x < out_size or y < out_size or x > (detector_sz - out_size) or y > (detector_sz - out_size): 
                pass
            # in bounds
            else:
                if uncertainty > 0:
                    mu, sigma = 0, uncertainty # mean and standard deviation

                    # generate random radius r with standard deviation sigma in mm
                    r = np.random.normal(mu, uncertainty/np.sqrt(2))
                    #r = np.random.poisson(uncertainty/np.sqrt(2))
                    #r = np.random.uniform(-1*uncertainty/np.sqrt(2),uncertainty/np.sqrt(2))
                    azimuth = np.random.uniform(0,2*np.pi)

                    newx = (r/10) * np.cos(azimuth) 
                    newy = (r/10) * np.sin(azimuth)

                    # make sure still in bounds
                    x = x + newx
                    y = y + newy
                    
                    if x < out_size or y < out_size or x > (detector_sz - out_size) or y > (detector_sz - out_size):
                        pass 
                    else:
                        xxes.append(x)
                        yxes.append(y)
                else:
                    xxes.append(x)
                    yxes.append(y)

        else:
            low_energy_electrons+=1
    #plt.hist(rr,bins=50)
    #plt.show()
    #plt.close()
    # ------------------------------- ------------------------------- ------------------------------- -------------------------------
    heatmap, xedges, yedges = np.histogram2d(xxes, yxes, bins=multiplier*nElements)

    fname_step = fname_save + '_raw.png'
    plot_step(heatmap,np.amax(heatmap),fname_step,label='# particles')

    # add noise here -- poisson noise
    # add isotropic low energy background(?)--one simulation 

    # first get the mask to use in 
    mask, decode = make_mosaic_MURA(nElements,boxdim,holes=False,generate_files=False)

    #fname_step = fname_save + '_mask1.png'
    #plot_step(mask,np.amax(mask),fname_step,label='# particles')

    decode = get_decoder_MURA(mask,nElements,holes_inv=False,check=check)
    decode = np.repeat(decode, multiplier, axis=1).repeat(multiplier, axis=0)

    #fname_step = fname_save + '_mask.png'
    #plot_step(decode,np.amax(decode),fname_step,label='# particles')

    # flip the heatmap over both axes bc point hole
    rawIm = np.fliplr(np.flipud(heatmap))

    # reflect bc correlation needs to equal convolution
    rawIm = np.fliplr(rawIm)

    # deconvolve
    result_image = fft_conv(rawIm,decode)
    fname_step = fname_save + '_dc.png'
    plot_step(result_image,np.amax(result_image),fname_step,label='# particles')
    #fig = plt.figure()
    #ax = plt.axes(projection='3d')

    #X = np.linspace(0, len(result_image), len(result_image));
    #Y = np.linspace(0, len(result_image), len(result_image));
    #X,Y = np.meshgrid(X,Y);

    #ax.plot_surface(X, Y, result_image, cmap='coolwarm')
    #plt.show()

    # take snr of only the noise floor - nope
    #if 'snr' in ff:
    #    nbins = multiplier*nElements
    #    std_result = []
    #    for nb in range(nbins):
    #        for nbi in range(nbins):
    #            # cut out middle section
    #            if nb / nbins > 0.4 and nb / nbins < 0.6:
    #                if nbi / nbins > 0.4 and nbi / nbins < 0.6:
    #                    pass
    #            else:
    #                std_result.append(result_image[nb,nbi])
    #  
    #    snr = np.amax(np.abs(result_image))/np.std(np.abs(np.array(std_result)))
    #else:
    #    snr = 0
    snr = np.amax(np.abs(result_image))/np.std(np.abs(result_image))

    # line plot of the diagonal -- must be flipped left and right but i have no clue why
    max_ind = np.where(result_image == np.amax(result_image))
    max_col = max_ind[1]
    if np.shape(max_col)[0] > 1:
        max_col = max_col[0]

    #resolution = plot_peak(np.fliplr(result_image)[:,int(max_col)],fname_save,condition='half_val')
    # for diagonal
    #resolution = plot_peak(np.diagonal(np.fliplr(result_image)),fname_save,condition='half_val')

    # for a moving diagonal
    #result_image = np.fliplr(result_image)
    #ind1 = (np.shape(result_image)[0]//2) - int(max_col)
    #result_image_cut = result_image[:,ind1:]
    #print(np.shape(result_image_cut))
    #tesolution = plot_peak(np.diagonal(result_image_cut),fname_save,condition='half_val')

    del xxes, yxes
    del result_image, heatmap, xedges,yedges, mask,decode,rawIm
    gc.collect()
    resolution = None
    
    return snr, resolution

dec_plt('data/hits.csv',0, 67, 0.66, 'testing_perp',87.78)