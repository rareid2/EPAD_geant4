import numpy as np
from plot_settings import *
import matplotlib.pyplot as plt 
import matplotlib as mpl
from matplotlib.patches import Ellipse, Circle 
import matplotlib.colors as colors
import seaborn as sns 
from os import listdir
from os.path import isfile, join

# ------------------------plotting function!-----------------------------
def plot_it(resolutions_list, labels, mt, colors):
    sns.set_palette("Paired")
    mask_distance = np.linspace(1,9,15) # cm
    mask_distance = mask_distance[1:] # cm
    positions_list = np.linspace(1,380,15)
    positions_list = positions_list[:-1]
    mask_distance = 2*np.rad2deg(np.arctan(positions_list/999.9))

    # fix figure size
    fig, ax = plt.subplots(figsize=(5,2.8))

    # Make the minor ticks and gridlines show.
    ax.grid(which='major', color='#DDDDDD', linewidth=0.5)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    plt.locator_params(axis="y", nbins=7)

    #ax2=ax.twinx()

    linestyles = ['dashed','dashed','dashed']
    ii=0
    for resolutions,label,color,linestyle in zip(resolutions_list,labels,colors,linestyles):
        ii+=1
        ax.plot(mask_distance,resolutions,label=label,linestyle=linestyle, marker='.',zorder=1,color=color)

        #ax.plot(mask_distance,resolutions[1],label=label,linestyle=linestyle, marker='.',zorder=1,color=color)
        #ax.fill_between(mask_distance, resolutions[0], resolutions[1], alpha=0.5,color=color)

        #if ii==2:
        #    ax.plot(mask_distance,resolutions[0],label=label,linestyle=linestyle, marker='.',zorder=1,color=color)
        #    ax.plot(mask_distance,resolutions[1],label=label,linestyle=linestyle, marker='.',zorder=1,color=color)
        #    ax.fill_between(mask_distance, resolutions[0], resolutions[1], alpha=0.5,color=color)
    
    #plt.xlabel('distance between mask and detector [cm]')
    ax.set_xlabel('geometric angle [deg]')

    #plt.xlabel('steps from center of mask')
    #plt.xlabel('position resolution [mm]')
    
    #plt.ylabel('fcfov [deg]')
    ax.set_ylabel('snr')
    #ax2.set_ylabel('ang res [deg]')

    #plt.ylabel('ang res [deg]')
    #plt.title(str(mt)+' um',color='#BFACC8')
    #plt.title('SNR vs noise figure',color='#BFACC8')
    #plt.legend()
    
    #plt.ylim([0,7])
    #plt.ylim([0,250])
    #plt.ylim([20,115])
    plt.ylim([0,80])
    #plt.ylim([0,7])
    #ax.set_ylim([5,50])
    #plt.xlim([1.25,9.25])
    #plt.xlim([-1,23])

    fname_save = 'results/parameter_sweeps/'+str(mt)+'_snr.png'

    fig.tight_layout(pad=0.1)
    plt.savefig(fname_save,dpi=300)
    plt.clf()

# import the data from the results files
all_res = []

res = []
data_file = 'results/parameter_sweeps/67_3300_6_0_snr.txt'
res1 = np.loadtxt(data_file)
#res1 = res1[:-2]

res.append(res1)

my_pallette = ['#e07a5f','#3d405b','#81b29a']
#"#9799CA"
colors = my_pallette[:3]

labels=['67']

# plot the results
plot_it(res, labels, '3300', [colors[1]])


# FOV RESULTS
#mask_distance = np.linspace(1,9,15) # cm
#2.219
#2.233
#2.238
# avg 2.23
#fov = np.rad2deg(2*np.arctan(2.23/np.array(mask_distance)))
#plot_it([fov], labels,'fov', ['#BFACC8'])


# uniform distribution snr
#snr1 = [236.55372128257451, 236.28630858014023, 234.42159622810397, 232.23835712451242, 228.05882238678097, 224.7131559673365, 218.6498333678197, 212.83702325573003, 205.63577878488604, 199.33598080067304, 194.90017029068414, 184.52101627432228, 178.15623166409807, 171.7910991904458, 163.18868604894482, 156.5615061853508, 149.39470063277398, 144.7566751086987, 137.047499355151, 131.35220457762594, 124.55294340554583, 119.20931719821552, 117.48370848641325, 109.09621409866155, 104.76050530170531, 100.06892561243151, 97.05649653559371, 94.01216409871459, 88.17110646192394, 87.12194458831475, 83.49049201922247, 80.50099208741867, 77.35731969208095, 75.66648569276407, 73.74287244889481, 72.8934264982409, 71.2056599361986, 67.83572338261544, 65.6174891692202, 64.95130185469986, 63.76355600122856, 61.22329825388917, 60.58924999220803, 58.19931626565736, 57.5145290056055, 55.78352459039376, 54.83578115530217, 52.82175036874279, 51.64477169299774, 50.81276313538434]
# normal distributed snr
#snr2 = [236.40361668833978, 234.48169234624206, 228.6276169686504, 222.8186964745949, 213.0770167635854, 204.9874704614915, 194.2426439159462, 184.39174360908876, 174.3956585999516, 164.0949319360858, 155.69575013464672, 144.31583332370627, 138.88691999865512, 129.08097752609706, 124.65090682379245, 117.24323582709745, 111.5600349850317, 107.09558858415855, 102.00039819572083, 95.38434403681646, 92.1444988121083, 89.99220434174721, 85.53612910315536, 82.05628053824067, 78.37619670775223, 76.01095787801077, 73.37277098579068, 68.65515848055426, 68.44861169086627, 65.95894830568523, 63.265515698757206, 61.09477265283589, 59.93599167825408, 58.892907140862384, 56.12624282173435, 55.40279246111978, 55.489101299169, 52.41124724082494, 50.32444442109429, 50.02966397981961, 49.65335895118743, 48.392371464990696, 46.84389454156122, 44.654138070657105, 42.62358882313906, 43.66073603639498, 41.85848661988628, 40.09621749783957, 39.987511865727946, 39.1126394826472]
