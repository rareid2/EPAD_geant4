import numpy as np 
import matplotlib.pyplot as plt 
import os
import math 

# HEADS UP TO CHANGE THE BEAM WIDTH!!!!!!!!!!!!!!


def align_vectors(a, b):
    #b = b / np.linalg.norm(b) # normalize a
    #a = a / np.linalg.norm(a) # normalize b
    v = np.cross(a, b)
    # s = np.linalg.norm(v)
    c = np.dot(a, b)

    v1, v2, v3 = v
    h = 1 / (1 + c)

    Vmat = np.array([[0, -v3, v2],
                  [v3, 0, -v1],
                  [-v2, v1, 0]])

    R = np.eye(3, dtype=np.float64) + Vmat + (Vmat.dot(Vmat) * h)
    return R

# provide n particles, list of positions in cm, rotation angles (theta and phi in deg) and energies in keV
def create_macro(path, n_particles, positions, rotations, energies, world_size):
    path_to_macrofile = os.path.join(path, 'auto_run_beams.mac')

    with open(path_to_macrofile, 'w') as f:
        f.write('/run/numberOfThreads 40 \n')
        f.write('/run/initialize \n')
        f.write('/control/verbose 0 \n')
        f.write('/run/verbose 0 \n')
        f.write('/event/verbose 0 \n')
        f.write('/tracking/verbose 0 \n')

        for ii, (pos, rot, ene) in enumerate(zip(positions, rotations, energies)):
            pos_string = str(pos[0]) + ' ' + str(pos[1]) + ' ' + str(pos[2])

            f.write('/gps/particle e- \n')
            f.write('/gps/pos/type Plane \n')
            f.write('/gps/pos/shape Circle \n')
            f.write('/gps/pos/centre ' + pos_string + ' cm \n')
            
            if rot !=0:

                detector_loc = np.array([0,0,world_size*.45])
                src = np.array(pos)
                norm_d = np.linalg.norm(detector_loc - src)
                normal = (detector_loc - src) / norm_d

                z_ax = np.array([0,0,1])

                xprime = np.cross(normal, z_ax)
                xprime_normd = xprime/np.linalg.norm(xprime)
                yprime = np.cross(normal,xprime_normd)
                yprime_normd = yprime/np.linalg.norm(yprime)

                #rot_m = align_vectors(z_ax, normal)
                #print(rot_m)

                rot1_string =  str(xprime_normd[0]) + ' ' + str(xprime_normd[1])  + ' ' + str(xprime_normd[2])
                rot2_string =  str(yprime_normd[0]) + ' ' + str(yprime_normd[1])  + ' ' + str(yprime_normd[2])

                #rot1_string =  '1 ' + str(round(np.cos(np.deg2rad(rot[0])),3)) + ' ' + str(-1*round(np.sin(np.deg2rad(rot[0])),3))
                #rot2_string =  str(round(np.cos(np.deg2rad(rot[1])),3)) + ' 1 ' + str(-1*round(np.sin(np.deg2rad(rot[1])),3))

                f.write('/gps/pos/rot1 ' + rot1_string +  ' \n')
                f.write('/gps/pos/rot2 ' + rot2_string +  ' \n')

            f.write('/gps/ang/type iso \n')
            f.write('/gps/ang/mintheta 0 deg \n')
            f.write('/gps/ang/maxtheta 0.23 deg \n')
            #f.write('/gps/ang/maxtheta 0.08 deg \n')
            f.write('/gps/ang/minphi 0 deg \n')
            f.write('/gps/ang/maxphi 360 deg \n')
            f.write('/gps/energy ' + str(ene) + ' keV \n')
            f.write('/run/beamOn ' + str(n_particles) + ' \n')

        #f.write('/vis/scene/endOfEventAction accumulate <0> \n')
    f.close()

def find_disp_pos(theta, z_disp):

    # find displaced postion needed to get an angular displacement
    x_disp = z_disp * np.tan(np.deg2rad(theta))

    return x_disp
    
"""
# first, how far away is the 'infinite' plane? -- to the aperture or the detector?
# let's just assume aperture

aperture_size = 6.1 # cm
half_resolution = 3.8/2 # deg 
min_dist_cm = (aperture_size/2)/np.tan(np.deg2rad(half_resolution))

# world size is 
world_size = 175 # cm
# detector
detector_placement = world_size*.45 - 1.5 # cm

# placement of the source needs to be ..
src_z = -1*(math.ceil(detector_placement))

if np.abs(src_z) > world_size/2:
    print('error world needs to be re-built to be atleast '+ str(src_z*2) + ' long')
# TESTING
# now some settings for the macro file
abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"

# displacement
abs_path = os.getcwd()
macro_path = os.path.join(abs_path, 'macros')
n_particles = 1000
theta = 7 # deg
src_z = -503.95
world_size = 1111
detector_placement = world_size*.45 # cm

disp = find_disp_pos(theta, np.abs(src_z) + detector_placement)
disp = round(disp,2)
positions = [[0,0,src_z],[disp,disp,src_z]]
rotations = [0,1]
energies = [500,500]

#positions = [[0,0,src_z]]
#rotations = [[0,0]]
#energies = [500]

create_macro(macro_path, n_particles, positions, rotations, energies, world_size)
"""