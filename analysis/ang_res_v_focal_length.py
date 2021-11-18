import numpy as np 
import matplotlib.pyplot as plt 
from plot_settings import * 
import os

s = 0.01 # cm
f = np.linspace(0.1,6,100)
theta_min = 2*s / f
theta_min = theta_min * 180 / np.pi

# plot
plt.plot(f,theta_min,'r--')
plt.xlabel('distance between mask and detector [cm]')
plt.ylabel('min angular resolution of the PSD [deg]')

fname = 'results/fall21_results/ang_res_v_f.png'
fpath = os.getcwd()
fig_name = os.path.join(fpath, fname)
plt.savefig(fig_name)
