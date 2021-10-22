#!/usr/bin/python3.5


import pandas as pd
import numpy as np

from scipy.stats import norm, skewnorm

# Extracts and returns actual inital particle source angles
from .fnc_findSourceAngle import findSourceAngle

# read in the detector hits and extract useful info
def get_det_hits(fname):
    # Read in raw hit data
    detector_hits = pd.read_csv(fname,
                               names=["det","x", "y", "z","energy"],
                               dtype={"det": np.int8, "x":np.float64,
                               "y": np.float64, "z":np.float64, "energy":np.float64},
                               delimiter=',',
                                on_bad_lines='skip',
                               engine='c')

    n_entries = len(detector_hits['det'])

    if len(detector_hits['det']) == 0:
        raise ValueError('No particles hits on either detector!')
    elif 2 not in detector_hits['det']:
        raise ValueError('No particles hit detector 2!')

    deltaX = np.zeros(n_entries, dtype=np.float64)
    deltaZ = np.zeros(n_entries, dtype=np.float64)

    array_counter = 0
    energies = []
    for count, el in enumerate(detector_hits['det']):
        # pandas series can throw a KeyError if character starts line
        # TODO: replace this with parse command that doesn't import keyerror throwing lines
        while True:
            try:
                pos1 = detector_hits['det'][count]
                pos2 = detector_hits['det'][count+1]

                detector_hits['x'][count]
                detector_hits['z'][count]

                detector_hits['x'][count+1]
                detector_hits['z'][count+1]

            except KeyError:
                count = count + 1
                if count == n_entries:
                    break
                continue
            break

        # Checks if first hit detector == 1 and second hit detector == 2
        if np.equal(pos1, 1) & np.equal(pos2, 2):
            deltaX[array_counter] = detector_hits['x'][count+1] - detector_hits['x'][count]
            deltaZ[array_counter] = detector_hits['z'][count+1] - detector_hits['z'][count]
            energies.append(detector_hits['energy'][count])

            # Successful pair, continues to next possible pair
            count = count + 2
            array_counter = array_counter + 1
        else:
            # Unsuccessful pair, continues
            count = count + 1

    # Copy of array with trailing zeros removed
    deltaX_rm = deltaX[:array_counter]
    deltaZ_rm = deltaZ[:array_counter]

    del deltaX
    del deltaZ

    return detector_hits, deltaX_rm, deltaZ_rm, energies

# -----------------------------------------------------
def calculateAnglePerParticle(gap_in_cm):
    # call read in function
    detector_hits, deltaX_rm, deltaZ_rm, energies = get_det_hits('./data/hits.csv')

    # Find angles in degrees
    theta = np.rad2deg(np.arctan2(deltaZ_rm, gap_in_cm))
    phi = np.rad2deg(np.arctan2(deltaX_rm, gap_in_cm))

    # Fit a standard normal distribution to data
    try:
        x_theta = np.linspace(min(theta), max(theta))
        mu_theta, std_theta = norm.fit(theta)
        p_theta = norm.pdf(x_theta, mu_theta, std_theta)

        x_phi = np.linspace(min(phi), max(phi))
        mu_phi, std_phi = norm.fit(phi)
        p_phi = norm.pdf(x_phi, mu_phi, std_phi)

    except:
        pass

    # Fit skew normal distribution to data
    #TODO: write a check for sig_p RuntimeError when np.sqrt(-#)
    alpha_t, loc_t, scale_t = skewnorm.fit(theta)
    alpha_p, loc_p, scale_p = skewnorm.fit(phi)

    delta_t = alpha_t/np.sqrt(1+alpha_t**2)
    delta_p = alpha_t/np.sqrt(1+alpha_p**2)

    mean_t = loc_t + scale_t*delta_t*np.sqrt(2/np.pi)
    mean_p = loc_p + scale_p*delta_p*np.sqrt(2/np.pi)

    p_test = scale_p**2 * (1 - 2*(delta_p**2)/np.pi)
    if np.equal(0, np.round(p_test, 2)):
        sig_p = None
    else:
        sig_p = np.sqrt(p_test)

    t_test = scale_t**2 * (1 - 2*(delta_t**2)/np.pi)
    if np.equal(0, np.round(t_test, 2)):
        sig_t = None
    else:
        sig_t = np.sqrt(t_test)

    theta_actual, phi_actual, numberOfParticles = findSourceAngle()

    with open('./data/results.txt', 'a') as f:
        f.write(str(numberOfParticles) +
        ',' + str(theta_actual) + ',' + str(phi_actual) +
        ',' + str(round(np.mean(theta), 4)) + ',' + str(round(np.std(theta), 4)) +
        ',' + str(round(np.mean(phi), 4)) + ',' + str(round(np.std(phi), 4)) +
        ',' + str(round(np.median(theta), 4)) + ',' + str(round(np.median(phi), 4)) +
        ',' + str(round(mu_theta, 4)) + ',' + str(round(std_theta, 4)) +
        ',' + str(round(mu_phi, 4)) + ',' + str(round(std_phi, 4)) +
        ',' + str(round(mean_t,4)) + ',' + str(round(sig_t,4)) +
        ',' + str(round(mean_p,4)) + ',' + str(round(sig_p,4)) + '\n')
    
    avg_KE = np.average(np.array(energies))
    return theta, theta_actual, avg_KE

# ---------------------------------------------------------------------------------------
