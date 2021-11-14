#!/usr/bin/python3.5
import pandas as pd
import numpy as np
from scipy.stats import norm, skewnorm
import os
# Extracts and returns actual inital particle source angles
from fnc_find_source_angle import findSourceAngle
from fnc_get_det_hits import getDetHits

# -----------------------------------------------------
def calculateAnglePerParticle(gap_in_cm):
    path_to_data = os.path.join(os.path.dirname(os.path.abspath("data")),'data/')
    path_to_hits_file = os.path.join(path_to_data, 'hits.csv')

    # call read in function
    detector_hits, deltaX_rm, deltaZ_rm, energies = getDetHits(path_to_hits_file)

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
