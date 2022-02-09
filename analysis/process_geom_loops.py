import numpy as np 
import os
import pandas as pd 
import scipy.stats
from scipy.signal import convolve2d as conv2
from scipy.fft import fft2,ifft
from skimage.transform import resize
from scipy.ndimage import zoom 
import scipy.signal

from deconvolve_and_plot import dec_plt
# from nearly the size of the detector to
#  max size!
# from nearly the size of the detector to max size!
mask_size = [90] # mm

# from 400 to 3300 um --- thicknesses to check
mask_thickness = np.linspace(400,3300,5)

# vary distance a bit
#mask_distance = np.linspace(0.5,9,10) # cm
mask_distance = [4.0]
# running mosaicked masks
mask_rank = [67]

# run the following point src differentiation
# need to re-range this (5 to 15?)
#thetas = np.linspace(5,0.5,10)
thetas = np.linspace(3.5,4.5,5)

positions_list = np.linspace(1,120,5)

for ms in mask_size:
    for mt in mask_thickness:
        for rr in mask_rank:
            snrs = []
            for md in mask_distance:

                rank = rr*2-1
                thickness = int(mt) # make sure integer value
                distance = md

                # start w snr
                fname = 'data/mosaic_sim/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_zero.csv'
                #print(fname)
                uncertainty = '2unc'
                nElements = rr
                boxdim = round(ms/rank,4) # in mm

                #ff = str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_snr_'+uncertainty
                #snr = dec_plt(fname,uncertainty,nElements,boxdim,ff,ms)                
                #snrs.append(snr)

                for pi,po in enumerate(positions_list):
                    for ti, theta in enumerate(thetas):
                        fname = 'data/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_' + str(ti) + '_'+ str(pi)+'.csv'
                        
                        print(fname)
                        ff = str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2))+'_'+ str(ms) + '_'+str(ti)+'_'+ str(pi)+ '_'+uncertainty
                        snr_2 = dec_plt(fname,uncertainty,nElements,boxdim,ff,ms)

                #print(snrs)
