import numpy as np 
import scipy.stats
import os
from fnc_calc_angle_per_particle import calculateAnglePerParticle

def findTheoreticalDist(det1_thickness_um, gap_in_cm, charge_nmbr, rest_mass_kg):
    # particle angle from the simulation results
    thetas, theta_act, avg_KE = calculateAnglePerParticle(gap_in_cm)

    # characteristic length for silicon     
    X0 = .09370 # in meters
    # nope this is wrong - use this 
    # from https://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/014.html

    # charge number for ....
    z = charge_nmbr
    # thickness
    x = det1_thickness_um * 10**(-6) # meters
    print('characteristic length is:', x/X0)

    # rest mass -- gotta change this!
    m0 = rest_mass_kg
    # velocity
    c = 299792458 # m/s
    # conversion to Joules
    eV2J = 1.60218e-19 # joules
    # energy -- this is technically relativistic KE
    E = avg_KE*10**(6) * eV2J # joules
    # gamma
    gamma = 1 + (E / (m0*c**2))
    # beta
    beta = np.sqrt(1-(1/gamma)**2)
    # betac
    beta_c = beta * c
    # relativistic momentum
    p = gamma * m0 * beta_c

    # term 1 and 2 for the distribution
    t1 = (13.6 * 10**6 * eV2J / (beta_c*p)) 
    t2 = z*np.sqrt(x/X0)*(1+0.038*np.log(x * z**2 /(beta**2 * X0)))

    # distribution settings
    standard_deviation = np.rad2deg(t1*t2)
    mean = 0

    data_plt = thetas - np.ones(len(thetas))*theta_act

    # get x and y values for the theoretical distribution
    x_values = np.arange(min(thetas),max(thetas),0.1)
    y_values = scipy.stats.norm(mean, standard_deviation)

    return x_values, y_values