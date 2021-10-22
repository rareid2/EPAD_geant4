import numpy as np 
import matplotlib.pyplot as plt 
import scipy.stats
import os
from fncs.fnc_calc_angle_per_particle import calculateAnglePerParticle
from plot_settings import *

# quick script to compare results to theoretical scattering MCS distribution (guassian)

def find_theoretical_dist(thickness):
    # particle angle from the simulation results
    thetas, theta_act, avg_KE = calculateAnglePerParticle(3.0)

    # characteristic length for silicon     
    X0 = .09370 # in meters
    # nope this is wrong - use this 
    # from https://pdg.lbl.gov/2009/AtomicNuclearProperties/HTML_PAGES/014.html

    # charge number for electrons
    z = 2
    # thickness
    x = thickness * 10**(-6) # meters
    print('characteristic length is:', x/X0)

    # rest mass
    m0 = 6.64424e-27 # in kg
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

    # plot it!
    plt.plot(x_values, y_values.pdf(x_values), '--',color=CB91_Violet,label='theoretical')
    plt.hist(data_plt,density=True, bins=100,color=CB91_Blue,label='simulation results')
    plt.legend()
    plt.margins(x=0)
    plt.margins(y=0)
    plt.xlabel('deg')
    plt.title(str(int(avg_KE))+'MeV alpha scattering through detector 1  \n char length = ' + str(round(x/X0, 4)))
    # save directory
    folder_save = '/home/rileyannereid/workspace/geant4/Geant4_electron_detector/results/fall21_results/reproducing_distribution'
    fname = os.path.join(folder_save, 'distribution_results_alpha_'+str(int(avg_KE))+'MeV_'+str(thickness)+'um.png')
    plt.savefig(fname)
    plt.close()

thickness = 100 #um
find_theoretical_dist(thickness)
