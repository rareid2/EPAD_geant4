import numpy as np 
import matplotlib.pyplot as plt 
from plot_settings import * 
import os
import math 

from run_particle_beams import create_macro, find_disp_pos

# --------------------- constants ----------------------------
# world size is  --- keep this the same
world_size = 175 # cm
detector_placement = world_size*.45 - 1.5 # cm

# placement of the source needs to be ..
src_z = -1*(math.ceil(detector_placement))

# now some settings for the macro file
abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"
n_particles = 1000000

# --------------------- --------------------- 
# ranges of parameters (excluding thickness for now)
mask_thickness = np.linspace(400,3300,3) #um
mask_distance = np.linspace(0.5,6,5) #cm
mask_rank = [61, 101, 139, 181]# 223]
thetas = [7,6,5,4,3,2]

# after this, run again with theta = 1 for 400um first 
# then run the SNR ones over night (SNR on y axis, everything else the same)
for rr in mask_rank:
    for md in mask_distance:
        rank = rr
        thickness = 1850
        distance = md

        mask_element_size = round(63/rank,4) # in mm

        mura_filename = "../src/mask_designs/" + str(rank) + "MURA_matrix.txt"

        # create the config file 
        f  = open("/home/rileyannereid/workspace/geant4/EPAD_geant4/src/config.txt", "w+")
        f.write(str(rank) + '\n')
        f.write(str(thickness) + '\n')
        f.write(str(distance) + '\n')
        f.write(str(mask_element_size) + '\n')
        f.write(str(mura_filename) + '\n')

        f.close()

        # next, generate the macro file -- question how many thetas to run? at what resolution???? try this w 5 of each first

        # start with 8??
        for theta in thetas:
            disp = find_disp_pos(theta, np.abs(src_z) + detector_placement)
            disp = round(disp,2)

            positions = [[0,0,src_z],[disp,disp,src_z]]
            rotations = [[0,0],[-1.5*theta,-1.5*theta]]
            energies = [500,500]

            # create the macro file
            create_macro(abs_path, n_particles, positions, rotations, energies) 

            cwd = os.getcwd()

            # go to build to run this
            os.chdir('build')

            cmd = './main ../macros/auto_run_beams.mac'

            # run the simulation
            os.system(cmd)

            # load the results and re-name them
            os.rename('../data/hits.csv', '../data/hits'+str(rank)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(theta)+'.csv')

            os.chdir(cwd)

        
        print('finished run ' + str(rr) + ' and ' + str(md))

        