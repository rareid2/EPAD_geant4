#!/usr/bin/python3.5

import subprocess
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tqdm
import sys
import os
from fnc_calc_angle_per_particle import calculateAnglePerParticle

abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"


def generateAutoRunFile(theta_in_deg, phi_in_deg, n_particles, energy_in_keV):

    theta = np.deg2rad(theta_in_deg)
    phi = np.deg2rad(phi_in_deg)

    # starting position is 5 cm back
    # linearly sample the x and z dimensions
    y_offset = -5 # cm
    z_pos = y_offset * np.tan(theta)
    x_pos = abs(y_offset) * np.tan(phi)

    position_string = str(x_pos) + ' ' + str(y_offset) + ' ' + str(z_pos)

    # Calculate initial momentum direction for particle
    y_dir = 1
    z_dir = y_dir * np.tan(theta)
    x_dir = y_dir * np.tan(phi)

    dir_string = str(x_dir) + ' ' + str(y_dir) + ' ' + str(z_dir)
    path_to_macros = os.path.join(abs_path,'macros/')
    path_to_macrofile = os.path.join(path_to_macros, 'auto_run_file.mac')

    with open(path_to_macrofile, 'w') as f:
        f.write('/run/initialize \n')
        f.write('/control/verbose 1 \n')
        f.write('/run/verbose 1 \n')
        f.write('/event/verbose 1 \n')
        f.write('/tracking/verbose 1 \n')

        f.write('/gps/particle e- \n')
        f.write('/gps/energy ' + str(energy_in_keV) + ' keV \n')
        f.write('/gps/position ' + position_string  + ' \n')
        f.write('/gps/direction ' + dir_string + ' \n')
        f.write('/gps/pos/type Point \n')

        f.write('/run/beamOn ' + str(n_particles) + ' \n')



def executeAutoRunFile():
    # use absolute paths
    path_to_build = os.path.join(abs_path,'build')
    path_to_main = os.path.join(path_to_build,'main')

    path_to_macros = os.path.join(abs_path,'macros/')
    path_to_macrofile = os.path.join(path_to_macros, 'auto_run_file.mac')

    bashCommand = path_to_main+ " " + path_to_macrofile
    process = subprocess.Popen(bashCommand.split(), stdout = subprocess.PIPE)
    output, error = process.communicate()

    if error is not None:
        exception("Error in simulation")

    if output is None:
        exception("Error in simulation: no output")


def cleanDataDirectory():
    path_to_data = os.path.join(abs_path,'data/')
    path_to_hits_file = os.path.join(path_to_data, 'hits.csv')
    path_to_init_file = os.path.join(path_to_data, 'init_pos.csv')

    if os.path.isfile(path_to_hits_file):
        bashCleanCommand = 'rm ' + path_to_hits_file
        process = subprocess.Popen(bashCleanCommand.split(), stdout=subprocess.PIPE)
    if os.path.isfile(path_to_init_file):
        bashCleanCommand = 'rm ' + path_to_init_file
        process = subprocess.Popen(bashCleanCommand.split(), stdout=subprocess.PIPE)
    else:
        return

def writeHeaderLine():
    path_to_data = os.path.join(os.path.dirname(os.path.abspath("data")),'data/')
    path_to_results_file = os.path.join(path_to_data, 'results.txt')

    with open(path_to_results_file, 'w') as f:
        f.write('Number_Particles,Theta_actual,Phi_actual,Theta_mean,Theta_std,Phi_mean,Phi_std,'+
                'Theta_median,Phi_median,Theta_normfit,T_s_nf,Phi_normfit,P_s_nf,Theta_snormfit,T_s_snf,'+
                'Phi_snormfit,P_s_snf\n')



# Modified range to iterate over floats (i.e. 10.5 degrees, etc.)
def frange(start, stop, step):
     i = start
     while i < stop:
         yield i
         i += step


def main():

    #####################################
    ####### Edit these parameters #######
    #####################################
    energy = 500 #keV

    angle_resolution = 5
    min_angle = -45
    max_angle = 45

    numberOfParticles = 10000

    #####################################
    ####### ^^^^^^^^^^^^^^^^^^^^^ #######
    #####################################

    # Overwrites whatever results.txt file was already in directory
    #writeHeaderLine()

    stds = []

    with tqdm.tqdm(total=100, unit_scale=True) as pbar:
        for ai, angle in enumerate(frange(min_angle, max_angle, angle_resolution)):
        
            # rename the hits file to generate more
            if ai > 0: 
                path_to_data = os.path.join(abs_path,'data/')
                path_to_hits_file = os.path.join(path_to_data, 'hits.csv')
                path_to_hits_file_rn = os.path.join(path_to_data, 'hits_'+str(int(last_ai))+'.csv')
                
                os.rename(path_to_hits_file,path_to_hits_file_rn)
            
            pbar.set_postfix(angle=angle, refresh=True)

            # Generates auto run file given the parameters we wish to simulate over
            generateAutoRunFile(theta_in_deg=angle,
                                phi_in_deg=0,
                                n_particles=numberOfParticles,
                                energy_in_keV=energy)

            # Removes any raw hit files that were already in data directory 
            #cleanDataDirectory()

            # Runs simulation with autogenerated run file, outputs raw hit results into data directory
            executeAutoRunFile()

            # update this for naming purposes
            last_ai = ai

            # Processes raw hit data into statistical estimates, appends to results.txt
            #theta, theta_actual, avg_KE = calculateAnglePerParticle(3.0)

            # Progress bar update
            #pbar.update((max_angle-min_angle)/angle_resolution)

            #stds.append(np.std(theta))
    
if __name__=='__main__':
    main()
