import pathlib
import glob

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.style.use('apj.mplstyle')

fig = plt.figure(figsize=(6.5, 6.5))
ax1  = fig.add_axes((0.1, 0.5, 0.33, 0.33))
ax2  = fig.add_axes((0.55, 0.5, 0.33, 0.33))
ax3  = fig.add_axes((0.1, 0.1, 0.33, 0.33))
ax4  = fig.add_axes((0.55, 0.1, 0.33, 0.33))
cax = fig.add_axes((0.09, 0.96, 0.75, 0.04))


axs = [ax1, ax2, ax3, ax4]


vmin = -4
vmax = -3


csvfiles = glob.glob('*.csv')
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
sm.set_array([])

norm_g = matplotlib.colors.Normalize(vmin=vmin + 0.5, vmax=vmax)
sm_g = plt.cm.ScalarMappable(cmap='Greys', norm=norm)
sm_g.set_array([])



for i, fname in enumerate(csvfiles):
    if i  == len(axs):
        break
    ax= axs[i]
    name = fname.split('.csv')[0]

    data = np.genfromtxt(fname, delimiter=',', skip_header=1)
    M = data[:,0]
    FeH = data[:,1]
    r = data[:,3]

    ff, mm = np.meshgrid(np.unique(FeH), np.unique(M))
    rr = np.zeros_like(ff)
    rr[:] = np.nan
    for i in range(len(M)):
        rr[(mm == M[i])*(ff == FeH[i])] = r[i]

#    mm, tt = np.meshgrid(np.unique(mesh), np.unique(time))
#    rr = np.zeros_like(tt)
#    rr[:] = np.nan
#    for i in range(len(mesh)):
#        rr[(mm == mesh[i])*(tt == time[i])] = r[i]


    ax.pcolormesh(mm, ff, np.log10(rr), vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)

    for i in range(len(M)):
        ax.text(M[i], FeH[i], '{:.2f}'.format(np.log10(r[i])), ha='center', va='center', size=9, color=sm_g.to_rgba(np.log10(r[i])))

    
    ax.invert_xaxis()
    ax.set_xticks((0.9, 1.1, 1.3, 1.5, 1.7))
    ax.set_yticks((-1.2, -0.8, -0.4, 0.0, 0.4))
ax3.set_xlabel('Mass')
ax3.set_ylabel('Fe/H')

cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap('viridis'), norm=norm, orientation='horizontal')
cbar.set_label('r')

fig.savefig('mesa_r_spread.png', dpi=600, bbox_inches='tight')
fig.savefig('mesa_r_spread.pdf', dpi=600, bbox_inches='tight')


