#!/usr/bin/python3.5

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import time

def plotResults(energy):
    d = pd.read_csv('./data/results.txt')


    fig = plt.figure()

    # Mean
    plt.plot(d.index, d['Theta_mean'], color='b', label='Theta Estimate')
    plt.fill_between(d.index, d['Theta_mean']-d['Theta_std'], d['Theta_mean']+d['Theta_std'],
                     color='blue', alpha=0.2, label=None)

    # Median
    #plt.plot(d.index, d['Theta_median'], color='g')

    # Skew normal fit
    #plt.plot(d.index, d['Theta_snormfit'], color='purple', label=None)
    plt.fill_between(d.index, d['Theta_snormfit']-d['T_s_snf'], d['Theta_snormfit']+d['T_s_snf'],
                     color='orange', alpha=0.5, label=r'1 $\sigma$ confidence')
    plt.fill_between(d.index, d['Theta_snormfit']-2*d['T_s_snf'], d['Theta_snormfit']+2*d['T_s_snf'],
                     color='orange', alpha=0.25, label=r'2 $\sigma$ confidence')
    # Actual
    plt.plot(d['Theta_actual'], linestyle='-.', color='r', label='Theta Actual')


    n_particles = d['Number_Particles'][0]
    plt.text(min(d.index), max(d.index)+5, 'E = %i keV, N = %i particles ' % (int(energy), int(n_particles)) , fontsize=10)
    plt.ylabel('Angle (degrees)')
    plt.grid()
    plt.legend()
    plt.show()

    #path = '/home/grant/Documents/Research/Marshall_Research/electron_detector/results/'
    #fileName = 'angle_est_%i_kev_' % int(energy) + str(time.strftime("%m%d_%H%M%S"))
    #f = path + fileName + '.pdf'

    #fig.savefig(f, dpi=fig.dpi)

    #return fileName
