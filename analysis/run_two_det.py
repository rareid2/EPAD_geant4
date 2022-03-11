import numpy as np 
import matplotlib.pyplot as plt 
import os
import math 
from run_particle_beams import create_macro, find_disp_pos
from fnc_get_det_hits import getDetHits
import re
# --------------------- constants ----------------------------
# world size is ---- keep this the same as defined in detector construction! lets just get bigger why not
world_size = 1111 # cm
world_sizeXY = 2000 # cm
d_gap = 3.0 # cm

# now some settings for the macro file
abs_path = os.getcwd()
macro_path = os.path.join(abs_path, 'macros/auto_run_beams.mac')

# run million particles
n_particles = 1000000

# set up geometry
detector_placement = 2*world_size*.45 # cm
md = 0
src_z = -1*(detector_placement)-md

thicknesses = np.linspace(20,200,10)
energies = np.logspace(2,4,10)
uncertaintyies = [0.5,3.0]

for uncertainty in uncertaintyies:
    # start loops
    sigma_ene = []
    for energy_in_keV in energies:
        sigmas = []
        for thickness in thicknesses:
            # create the config file 
            f  = open("src/config.txt", "w+")
            # thickness in um
            f.write(str(thickness) + '\n')
            f.close()

            position_string = '0 0 ' + str(src_z)
            dir_string = '0 0 1'

            """
            # create macro file!
            with open(macro_path, 'w') as f:
                f.write('/run/numberOfThreads 40 \n')
                f.write('/run/initialize \n')
                f.write('/control/verbose 0 \n')
                f.write('/run/verbose 0 \n')
                f.write('/event/verbose 0 \n')
                f.write('/tracking/verbose 0 \n')

                f.write('/gps/particle e- \n')
                f.write('/gps/energy ' + str(energy_in_keV) + ' keV \n')
                f.write('/gps/position ' + position_string  + ' \n')
                f.write('/gps/direction ' + dir_string + ' \n')
                f.write('/gps/pos/type Point \n')

                f.write('/run/beamOn ' + str(n_particles) + ' \n')
                f.close()

            # execute!
            
            # go to build to run this
            os.chdir('build')
            cmd = './main ../macros/auto_run_beams.mac'
            # run the simulation
            os.system(cmd)
            # load the results and re-name them
            fname = 'data/hits_'+str(thickness)+'_'+str(energy_in_keV)+'.csv'
            os.rename('../data/hits.csv', '../'+fname)
            os.chdir(abs_path)
            """
            
            # process
            fname = 'data/hits_'+str(thickness)+'_'+str(energy_in_keV)+'.csv'
            print(fname)
            detector_hits, deltaX_rm, deltaZ_rm, energies_data = getDetHits(fname,uncertainty)
            thetas = np.arctan(deltaZ_rm / d_gap)
            one_sigma = np.std(np.rad2deg(thetas))

            # save
            sigmas.append(one_sigma)
            #print('finished one round of energies')
        sigma_ene.append(sigmas)
        print(sigmas)
    print(sigma_ene)