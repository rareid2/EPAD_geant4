import numpy as np 
import matplotlib.pyplot as plt 
import os
import math 
from run_particle_beams import create_macro, find_disp_pos
from deconvolve_and_plot import dec_plt
import re
# --------------------- constants ----------------------------
# world size is ---- keep this the same as defined in detector construction! lets just get bigger why not
world_size = 1111 # cm
world_sizeXY = 2000 # cm

# now some settings for the macro file
abs_path = os.getcwd()
macro_path = os.path.join(abs_path, 'macros')

# run million particles
n_particles = 1000000

detector_placement = world_size*.45 # cm

uncertainty = 0
unc = uncertainty
# --------------------- --------------------- 
# ranges of parameters!!

# from nearly the size of the detector to max size!
ms = 87.78 # mm

# from 400 to 3300 um --- thicknesses to check
mask_thickness = np.linspace(400,3300,3)
mask_thickness = [1850]

# vary distance a bit
mask_distance = np.linspace(1,9,15) # cm
mask_distance = mask_distance[1:]
mask_distance = mask_distance[:6]
mask_distance = [mask_distance[0]]

# running mosaicked masks
mask_rank = [67]

positions_list = np.linspace(1170,1480,10)
#positions_list = positions_list[:6]
#thetas = np.linspace(1.0,2.3,1)

#positions_list = positions_list[-2:]
# loop through everything
for ii, mt in enumerate(mask_thickness):
    for rr in mask_rank:
        #if rr == 31:
        #    ms = 26.84
        #else:
        #    ms = 26.62
        snrs = []
        res_dists = []

        for md in mask_distance:  
            #if md < 3.5:
            #    thetas = np.linspace(2,6,18)
            #else:
            #    thetas = np.linspace(0.75,2.25,15)
                
            # placement of the source needs to be ..
            src_z = -1*(detector_placement)-md         
            
            # get pixel size
            rank = rr*2-1
            thickness = int(mt) # make sure integer value
            distance = md
            nElements = rr
            boxdim = round(ms/rank,4) # in mm

            # find the pixel size
            mask_element_size = round(ms/rank,4) # in mm

            # load in the filename for the mask
            mura_filename = "../src/mask_designs/" + str(rr) + "mosaicMURA_matrix_"+str(ms)+".txt"
            
            # create the config file 
            f  = open("src/config.txt", "w+")

            # give it the total number of elements!
            f.write(str(rank) + '\n')

            # thickness in um
            f.write(str(thickness) + '\n')

            # distance in cm
            f.write(str(distance) + '\n')

            # size of element in mm 
            f.write(str(mask_element_size) + '\n')

            # filename for mask
            f.write(str(mura_filename) + '\n')

            f.close()

            for pi,po in enumerate(positions_list):
                #po = 0
                po = round(po,3)
                
                # first run the SNR
                positions = [[po,0,src_z]]
                rotations = [1]
                energies = [500]
                
                # create the macro file
                create_macro(macro_path, n_particles, positions, rotations, energies, world_size) 
                
                cwd = os.getcwd()
                # go to build to run this
                os.chdir('build')
                cmd = './main ../macros/auto_run_beams.mac'
                # run the simulation
                os.system(cmd)
                # load the results and re-name them
                os.rename('../data/hits.csv', '../data/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+ '_' + str(po)+ '_zero.csv')
                #os.rename('../data/hits.csv', '../data/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+'_zero.csv')
                
                os.chdir(cwd)
                
                ff = str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+  '_' + str(po)+ '_' +str(unc)
                fname = 'data/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(round(distance,2))+ '_' + str(po)+  '_zero.csv'
                snr, resolution2 = dec_plt(fname,uncertainty,nElements,boxdim,ff,ms)                
                snrs.append(snr)

                """
                for ti,theta in enumerate(thetas):

                    disp = find_disp_pos(theta, np.abs(src_z) + (detector_placement-md))
                    disp = round(disp,2)

                    positions = [[po,0,src_z],[po+disp,0,src_z]]
                    #positions = [[0,0,src_z],[disp,0,src_z]]

                    rotations = [1,1]
                    energies = [500,500]

                    # create the macro file 
                    
                    create_macro(macro_path, n_particles, positions, rotations, energies, world_size) 
                    
                    cwd = os.getcwd()
                    # go to build to run this
                    os.chdir('build')
                    cmd = './main ../macros/auto_run_beams.mac'
                    # run the simulation
                    os.system(cmd)
                    # load the results and re-name them
                    os.rename('../data/hits.csv', '../data/hits_'+str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2))+'_' + str(theta)+'_' + str(po)+'.csv')
                    #os.rename('../data/hits.csv', '../data/hits_'+str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2))+'_'+ str(ms) + '_' + str(theta)+'.csv')
                    os.chdir(cwd)
                    
                    # process it
                    fname = 'data/hits_'+str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2))+'_'+str(theta)+'_' + str(po) + '.csv'
                    ff = str(rr)+'_'+str(thickness)+'_'+ str(round(distance,2)) + '_'+str(po) +'_'+str(theta)+'_' +str(unc)
                    snr_2, resolution = dec_plt(fname,uncertainty,nElements,boxdim,ff,ms)

                    if resolution==True:
                        res_dists.append(thetas[ti])
                        break
                    else:
                        if thetas[ti] == max(thetas):
                            res_dists.append(max(thetas)+1)
                            continue
                        else:
                            continue
                """

        # save the results
        #a_file = open('results/parameter_sweeps/'+str(rr)+'_'+str(thickness) + '_' + str(md) + '_' +str(unc) + '_snr.txt', "w")
        #np.savetxt(a_file,snrs)
        #a_file.close()

        #a_file = open('results/parameter_sweeps/'+str(rr)+'_'+str(thickness)+'_' + str(md)+ '_' +str(unc) +'_res.txt', "w")
        #np.savetxt(a_file,res_dists)
        #a_file.close()