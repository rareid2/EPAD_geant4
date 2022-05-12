import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import AxesGrid
import seaborn as sns
from matplotlib import ticker, cm
# set plot
fig, axes = plt.subplots(3,3, figsize=(8,6),constrained_layout=True,sharex=True, sharey=True)

L = 62 #mm
Rl = 9.92 #kOhm

# vary temperature from 10 to 300
Ts = np.linspace(10,300,1000)
# vary energy deposited from 100 to 1000
E_deps = np.linspace(50,1000,1000) #keV

# create grid space
E_dep, T = np.meshgrid(E_deps, Ts)

# vary equivalent resistance from 1kOhm to 1MOhm?
Reqs = [1,100,1000]

# vary tau from 1us to 10us
taus = [1e-6,5e-6,10e-6]

cmap_reversed = cm.get_cmap('jet')

for i,Req in enumerate(Reqs):
    for j,tau in enumerate(taus):
        delat_L = 14.8 * (L/E_dep) * np.sqrt( (T * tau) * (1 + (Req/Rl) )/ (Rl/2))
        cs = axes[i,j].pcolormesh(E_dep, T, delat_L, cmap=cmap_reversed,vmin=0,vmax=5)


fig.supylabel('Temperature [K]')
fig.supxlabel('Energy Deposited [keV]')
cbar = fig.colorbar(cs, ax=axes[:, 2],fraction=0.6,label='position uncertainty [mm]')
plt.savefig('testing.png')