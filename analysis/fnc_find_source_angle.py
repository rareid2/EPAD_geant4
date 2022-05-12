#!/usr/bin/python3.5
import pandas as pd
import numpy as np

def findSourceAngle():

    initial_params = pd.read_csv('./data/init_pos.csv',
                                 names=["x","y","z","momX","momY","momZ"],
                                 dtype=np.float64)

    # Only samples first particle since all have same direction from a point source
    theta_actual = round(np.rad2deg(np.arctan2(initial_params['momZ'][0], initial_params['momY'][0])), 4)
    phi_actual = round(np.rad2deg(np.arctan2(initial_params['momX'][0], initial_params['momY'][0])), 4)

    numberOfParticles = int(initial_params.size/6)

    return theta_actual, phi_actual, numberOfParticles
