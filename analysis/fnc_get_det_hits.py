import pandas as pd
import numpy as np

# read in the detector hits and extract useful info
def getDetHits(fname):
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