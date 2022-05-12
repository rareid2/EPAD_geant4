import numpy as np 
import scipy.stats
import os
from fnc_calc_angle_per_particle import calculateAnglePerParticle
from fnc_get_det_hits import getDetHits

def findTheoreticalDist(det1_thickness, gap_in_cm, charge_nmbr, rest_mass_ME,avg_KE,data_path=None):
    """
    if data_path == None:
        path_to_data = os.path.join(os.path.dirname(os.path.abspath("data")),'data/')
        data_path = os.path.join(path_to_data, 'hits.csv')

    # particle angle from the simulation results
    detector_hits, deltaX_rm, deltaZ_rm, energies = getDetHits(data_path)
    avg_KE = np.average(energies)
    thetas = np.rad2deg(np.arctan2(deltaX_rm, gap_in_cm))
    """
    #avg_KE = 1000

    X0 = 352.7597513557388*10**3 # in um
    #X0 = 37*10**3 # in um

    # charge number for ....
    z = charge_nmbr
    # thickness
    x = det1_thickness

    #print('characteristic length is:', x/X0)
    if x/X0 > 100:
        print(x,avg_KE)

    # updated beta_cp (from src code)
    E = avg_KE * 0.001 # convert to ME
    invbeta_cp = (E + rest_mass_ME)/(E**2 + 2*rest_mass_ME*E) 

    # term 1 and 2 for the distribution -- urban 2006 eqn 28
    Z_si = 4
    fz = 1 - (0.24/(Z_si*(Z_si + 1)) )
    t1 = 13.6 * invbeta_cp # in ME
    t2 = z*np.sqrt(x/X0)*np.sqrt(1+0.105*np.log(x/X0) + 0.0035*(np.log(x/X0))**2)*fz

    # distribution settings
    standard_deviation = np.rad2deg(t1*t2)
    mean = 0

    # get x and y values for the theoretical distribution
    #x_values = np.arange(min(thetas),max(thetas),0.1)
    x_values = np.arange(-60,60,0.1)
    
    y_values = scipy.stats.norm(mean, standard_deviation)
    #print(standard_deviation)

    return x_values, y_values, standard_deviation