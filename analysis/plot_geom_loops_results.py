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
def plot_it(resolutions_list, labels, mt, color):
    sns.set_palette("Paired")
    mask_distance = np.linspace(0.5,6,5) #cm
    
    # fix figure size
    fig, ax = plt.subplots(figsize=(5,2.8))

    for resolutions,label in zip(resolutions_list,labels):
        plt.plot(mask_distance,resolutions,label=label,linestyle='dashed', marker='s')
    
    plt.legend(loc='upper right')
    plt.xlabel('distance between mask and detector [cm]')
    plt.ylabel('angular resolution [deg]')
    plt.title(str(mt)+'um',color=color)
    plt.ylim([10,75])

    abs_path = "/home/rileyannereid/workspace/geant4/EPAD_geant4"
    fname_save = 'results/fall21_results/ptsrcdiff/'+str(mt)+'_snr_unc.png'

    fig.tight_layout(pad=0.1)
    plt.savefig(fname_save)
    plt.clf()

# parameter ranges
mask_thickness = np.linspace(400,3300,3) #um
mask_distance = np.linspace(0.5,6,5) #cm
mask_rank = [61,101,139,181]
labels = ['61','101','139','181']
colors = ['#FFDC5E', '#DC493A', '#CB9CF2']

# note that only integer values possible
resolutions_list = []

# for 400um
#res_61 = [8,6,4,3,3]
#res_101 = [8,4,3,1,2]
#res_139 = [8,4,3,2,3]
#res_181 = [8,3,2,3,4]

# for 400um w uncertainty
#res_61 = [8,8,5,4,3]
#res_101 = [8,8,5,4,3]
#res_139 = [8,8,5,4,3]
#res_181 = [8,6,4,4,3]

# for 1850um
#res_61 = [8,6,4,3,3]
#res_101 = [8,4,3,3,2]
#res_139 = [8,4,3,3,3]
#res_181 = [8,3,2,3,3]

# for 1850um w uncertainty
#res_61 = [8,8,5,4,3]
#res_101 = [8,8,5,4,3]
#res_139 = [8,8,5,4,3]
#res_181 = [8,7,5,4,3]

# for 3300um 
#res_61 = [8,7,4,3,3]
#res_101 = [8,5,3,2,2]
#res_139 = [8,4,3,3,3]
#res_181 = [8,3,3,3,3]

# for 3300um w uncertainty
#res_61 = [8,8,6,4,3]
#res_101 = [8,8,6,4,3]
#res_139 = [8,8,5,4,3]
#res_181 = [8,7,5,5,3]

# SNR
# 400um
res_61 = [54.86918689072489, 52.51162116079835, 45.381033372515915, 34.03872948181812, 24.66703898631042]
res_101 = [74.76376626699493, 61.99754121536994, 37.800493676521704, 21.48060378345827, 15.399927374040669]
res_139 = [73.61874113304235, 49.56877946068146, 26.61356018364926, 17.14249246168529, 16.742452520098706]
res_181 = [63.03732030456078, 39.38148860450789, 21.13707297961793, 12.245596317080794, 11.260043012758356]

# 1850um
res_61 = [54.925413971890215, 52.70438970612835, 45.57567550656722, 34.684714722889495, 25.287101671758986]
res_101 = [74.85153758835942, 63.276931328669804, 39.041547281720796, 22.054314027788887, 15.217523510901648]
res_139 = [74.72716663132432, 51.95822513023131, 27.468386087240262, 17.866358476995774, 16.65592352822767]
res_181 = [65.00638160245897, 42.046411770870144, 21.994475740030584, 12.38721964438829, 11.032613989677284]


# 3300um
res_61 = [54.80934678956001, 53.004660477157195, 46.11269023871586, 35.29992456557015, 25.537422085005566]
res_101 = [75.26808814558342, 64.26244959966719, 40.55008048662667, 23.097999929535757, 14.84991942056127]
res_139 = [75.92172412437404, 54.06934313982156, 28.82584889004756, 17.818196873033976, 16.796538109387345]
res_181 = [66.13270296955385, 44.6048049824886, 23.3297086278307, 12.03525047484963, 11.136965452465635]


# with uncertainty
# 400
res_61 = [28.46439872728806, 32.226020279191616, 34.53322312703455, 34.568898718915406, 34.82456726859576]
res_101 = [17.60005097690467, 19.564625462930845, 22.857753227520707, 22.427730233721196, 23.247480086189306]
res_139 = [14.0139240405507, 14.937421901038388, 16.20987585643321, 14.968469108792284, 15.699913570798126]
res_181 = [10.68641774848388, 12.709144936490388, 13.447124999263707, 13.22213453668203, 12.636234230492326]

# 1850
res_61 = [28.53797599647956, 32.45460469693967, 34.1688905159045, 35.46286152404505, 33.81238077261387]
res_101 = [17.867286947113875, 19.446914865625338, 21.59274692459884, 21.57987465660968, 21.093712184783403]
res_139 = [15.175910197275405, 15.367224454650616, 16.465194922401416, 17.533795022566025, 16.489064610154436]
res_181 = [10.722457021440963, 11.922995636381268, 12.131644830437425, 12.906542936654592, 12.871437742547998]

# 3300
res_61 = [28.492328411585678, 32.74742175758, 34.74118131715422, 34.53361682261932, 34.010666389087966]
res_101 = [20.16846867120859, 21.240231486917207, 22.749280756314132, 21.945145304624905, 22.49972059910511]
res_139 = [15.77608337088403, 16.352284403938654, 15.55804531513328, 16.54245414828009, 16.69873272283704]
res_181 = [12.160887188496043, 12.258723484625182, 14.310138710244141, 13.016145032403555, 13.186880988830563]


resolutions_list.append(res_61)
resolutions_list.append(res_101)
resolutions_list.append(res_139)
resolutions_list.append(res_181)

mask_thickness = np.linspace(400,3300,3) #um
mask_distance = np.linspace(0.5,6,5) #cm
mask_rank = [61,101,139,181]

plot_it(resolutions_list, labels, 3300, colors[2])


"""
# start with theoretical first
for mi,mt in enumerate(mask_thickness):
    th_resolutions = []
    labels = []
    for mr in mask_rank:
        pixel_size = round(63/mr,4)
        th_res = 2*np.rad2deg(np.arctan((pixel_size/2)/(np.array(mask_distance - (mt*1e-4))*10)))

        th_resolutions.append(th_res)
        labels.append(str(mr))

    plot_it(th_resolutions, labels, round(mt), colors[mi])
"""