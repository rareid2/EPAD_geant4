#!/usr/bin/python3.5

from astroML.stats import fit_bivariate_normal
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

detector1_hits = pd.read_csv("./data/hits_det1.csv",
                             names=["x", "y", "z", "mom_x", "mom_y", "mom_z"],
                             dtype=np.float64)

detector2_hits = pd.read_csv("./data/hits_det2.csv",
                             names=["x", "y", "z", "mom_x", "mom_y", "mom_z"],
                             dtype=np.float64)

# Hit data from detector 1
X1 = detector1_hits["x"]
Y1 = detector1_hits["z"]

# Hit data from detector 2
X2 = detector2_hits["x"]
Y2 = detector2_hits["z"]

###########################################
### astroML's robust fit for detector 1 ###
###########################################

(mu_nr1, sigma1_nr1,
     sigma2_nr1, alpha_nr1) = fit_bivariate_normal(X1, Y1, robust=True)

fig = plt.figure()
ax = fig.add_subplot(1,2,1)

# Scatter plot of hit data
ax.scatter(X1, Y1, s=2, lw=0, c='k', alpha=0.25)

E_2sigma = Ellipse(mu_nr1, sigma1_nr1 * 2, sigma2_nr1 * 2,
            (alpha_nr1 * 180. / np.pi), ec='b', fc='none', linestyle='dotted')

ax.add_patch(E_2sigma)

E_3sigma = Ellipse(mu_nr1, sigma1_nr1 * 3, sigma2_nr1 * 3,
            (alpha_nr1 * 180. / np.pi), ec='r', fc='none', linestyle='dotted')

ax.add_patch(E_3sigma)

ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_xlabel("x [cm]")
ax.set_ylabel("z [cm]")
ax.text(0.04, 0.96, "Detector 1",
            ha='left', va='top', transform=ax.transAxes)

ax.plot([0], [0], '-b', label='2 sigma')
ax.plot([0], [0], '-r', label='3 sigma')
ax.legend(loc='lower right')


###########################################
### astroML's robust fit for detector 2 ###
###########################################

(mu_nr2, sigma1_nr2,
     sigma2_nr2, alpha_nr2) = fit_bivariate_normal(X2, Y2, robust=True)

ax = fig.add_subplot(1,2,2)

# Scatter plot of hit data
ax.scatter(X1, Y1, s=2, lw=0, c='k', alpha=0.25)

E_2sigma = Ellipse(mu_nr2, sigma1_nr2 * 2, sigma2_nr2 * 2,
            (alpha_nr2 * 180. / np.pi), ec='b', fc='none', linestyle='dotted')

ax.add_patch(E_2sigma)

E_3sigma = Ellipse(mu_nr2, sigma1_nr2 * 3, sigma2_nr2 * 3,
            (alpha_nr2 * 180. / np.pi), ec='r', fc='none', linestyle='dotted')

ax.add_patch(E_3sigma)

ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_xlabel("x [cm]")
ax.set_ylabel("z [cm]")
ax.text(0.04, 0.96, "Detector 2",
            ha='left', va='top', transform=ax.transAxes)

ax.plot([0], [0], '-b', label='2 sigma')
ax.plot([0], [0], '-r', label='3 sigma')
ax.legend(loc='lower right')

plt.show()
