import numpy as np 
import matplotlib.pyplot as plt 
from plot_settings import * 
import os
import math 
from run_particle_beams import create_macro, find_disp_pos

# --------------------- constants ----------------------------
# world size is ---- keep this the same as defined in detector construction! lets just get bigger why not
world_size = 250 # cm

# now some settings for the macro file
abs_path = os.getcwd()
macro_path = os.path.join(abs_path, 'macros')
# run million particles
n_particles = 1000000

detector_placement = world_size*.45 # cm
# placement of the source needs to be ..
src_z = -1*(math.ceil(detector_placement))

# --------------------- --------------------- 
# ranges of parameters!!

# from nearly the size of the detector to max size!
mask_size = [73,82,90] # mm

# from 400 to 3300 um --- thicknesses to check
mask_thickness = [400,1850,3300] # um

# vary distance a bit
mask_distance = np.linspace(0.5,9,10) # cm

# running mosaicked masks
mask_rank = [37, 67, 97, 127]

# run the following point src differentiation
thetas = np.linspace(6,0.5,10)

# loop through everything
for ms in mask_size:
    for mt in mask_thickness:
        for rr in mask_rank:
            for md in mask_distance:
                
                # get pixel size
                rank = rr*2-1
                thickness = int(mt) # make sure integer value
                distance = md
                size = ms

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

                # first run the SNR
                positions = [[0,0,src_z]]
                rotations = [[0,0]]
                energies = [500]

                # create the macro file
                create_macro(macro_path, n_particles, positions, rotations, energies) 

                cwd = os.getcwd()
                # go to build to run this
                os.chdir('build')
                cmd = './main ../macros/auto_run_beams.mac'
                # run the simulation
                os.system(cmd)
                # load the results and re-name them
                os.rename('../data/hits.csv', '../data/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(ms) + '_zero.csv')
                os.chdir(cwd)


                for ti,theta in enumerate(thetas):
                    disp = find_disp_pos(theta, np.abs(src_z) + detector_placement)
                    disp = round(disp,2)

                    positions = [[0,0,src_z],[disp,disp,src_z]]
                    rotations = [[0,0],[-1.5*theta,-1.5*theta]]
                    energies = [500,500]

                    # create the macro file
                    create_macro(macro_path, n_particles, positions, rotations, energies) 

                    cwd = os.getcwd()
                    # go to build to run this
                    os.chdir('build')
                    cmd = './main ../macros/auto_run_beams.mac'
                    # run the simulation
                    os.system(cmd)
                    # load the results and re-name them
                    os.rename('../data/hits.csv', '../data/hits_'+str(rr)+'_'+str(thickness)+'_'+  str(distance)+'_'+ str(ms) + '_' + str(ti)+'.csv')

                    os.chdir(cwd)