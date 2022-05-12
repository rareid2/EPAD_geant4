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

mask_size = [90] # mm
ms = mask_size[0]

# from 400 to 3300 um --- thicknesses to check
mask_thickness = [400,1850,3300]

# vary distance a bit
mask_distance = np.linspace(1,9,10) # cm
mask_distance = [5.44]
# running mosaicked masks
mask_rank = [37,67,97]

# run the following point src differentiation
thetas = np.linspace(0.75,5.5,10)

for mt in mask_thickness:
    for rr in mask_rank:
        snrs = []
        res_dists = []
        for md in mask_distance:

            rank = rr*2-1
            thickness = int(mt) # make sure integer value
            distance = md

            # start w snr
            fname = 'data/mosaic_sim3/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_zero.csv'
            print(fname)
            
            uncertainty = 0
            nElements = rr
            boxdim = round(ms/rank,4) # in mm

            ff = str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_snr_'+str(uncertainty)
            snr, resolution2 = dec_plt(fname,uncertainty,nElements,boxdim,ff,ms)                
            snrs.append(snr)

            #for pi,po in enumerate(positions_list):
            resolutions = []
            for ti, theta in enumerate(thetas):
                #fname = 'data/mosaic_sim3/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_' + str(ti) + '_'+ str(pi)+'.csv'
                fname = 'data/mosaic_sim3/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_'+ str(ms) + '_' + str(ti) + '.csv'
                
                print(fname)
                #ff = str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2))+'_'+ str(ms) + '_'+str(ti)+'_'+ str(pi)+ '_'+uncertainty
                ff = str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2))+'_'+ str(ms) + '_'+str(ti)+'_'+str(uncertainty)
                
                snr_2, resolution = dec_plt(fname,uncertainty,nElements,boxdim,ff,ms)

                if resolution==True:
                    resolutions.append(thetas[ti])
                else:
                    resolutions.append(6)

            res_dist = np.min(np.array(resolutions))
            res_dists.append(res_dist)
        print(snrs)
        print(res_dist)
