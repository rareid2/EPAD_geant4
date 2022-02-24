import pandas as pd
import numpy as np

# read in the detector hits and extract useful info
def getDet1Hits(fname):
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

    posX = []
    posY = []
    energies = []
    for count, el in enumerate(detector_hits['det']):
        # only get hits on the first detector
        if el == 1:
            xpos = detector_hits['x'][count]
            zpos = detector_hits['y'][count]
            energy_kev = detector_hits['energy'][count]

            # save
            posX.append(xpos)
            posY.append(zpos)
            energies.append(energy_kev)

    # remove(?) helps w memory(?)
    del detector_hits
    
    return posX, posY, energies
