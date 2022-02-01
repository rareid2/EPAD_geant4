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

# mask size will vary
mask_size = [30,60,90] # mm
mask_size = [84]
# from 400 to 3300 um --- thicknesses to check
mask_thickness = [400,1850,3300] # um
mask_thickness = [400]
# vary distance a bit
mask_distance = np.linspace(0.5,9,10) # cm
mask_distance = [2.0]
# running mosaicked masks
mask_rank = [37, 67, 97]
mask_rank = [11]
# run the following point src differentiation
thetas = np.linspace(6,0.5,10)

for ms in mask_size:
    for mt in mask_thickness:
        for rr in mask_rank:
            for md in mask_distance:

                rank = rr*2-1
                thickness = int(mt) # make sure integer value
                distance = md

                # start w snr 
                fname = 'data/save_hits/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(ms) + '_zero.csv'
                uncertainty = False
                nElements = rr
                boxdim = round(ms/rank,4) # in mm

                ff = str(rr)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(ms) + '_snr'
                snr = dec_plt(fname,uncertainty,nElements,boxdim,ff)

                #for ti, theta in enumerate(thetas):
                #    fname = 'data/save_hits/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(ms) + '_' + str(ti)+'.csv'
                #    ff = str(rr)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(ms) + '_'+str(ti)
                #    snr = dec_plt(fname,uncertainty,nElements,boxdim,ff)
