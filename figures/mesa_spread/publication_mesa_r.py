import pathlib
import glob

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.style.use('apj.mplstyle')

fig = plt.figure(figsize=(6.5, 6.5))
ax1  = fig.add_axes((0.1, 0.55, 0.35, 0.35))
ax2  = fig.add_axes((0.55, 0.55, 0.35, 0.35))
ax3  = fig.add_axes((0.1, 0.1, 0.35, 0.35))
ax4  = fig.add_axes((0.55, 0.1, 0.35, 0.35))
cax = fig.add_axes((0.1, 0.96, 0.80, 0.04))

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

for fname in csvfiles:
    if 'Brown' in fname:
        i = 0
        tag = 'Brown'
    elif 'coeff1.0e-01' in fname:
        i = 1
        tag = r'Kipp, $\alpha_{\rm{th}} = 0.1$'
    elif 'coeff2.0e+00' in fname:
        i = 2
        tag = r'Kipp, $\alpha_{\rm{th}} = 2$'
    elif 'coeff7.0e+02' in fname:
        i = 3
        tag = r'Kipp, $\alpha_{\rm{th}} = 700$'

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

    bbox_props = dict(ec="k", alpha=1, fc="w", linewidth=0.8, pad=2)
#    bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.75)
    ax.text(0.0128, 1.0125, tag, ha="left", va="bottom", size=11, transform=ax.transAxes, bbox=bbox_props)

#    mm, tt = np.meshgrid(np.unique(mesh), np.unique(time))
#    rr = np.zeros_like(tt)
#    rr[:] = np.nan
#    for i in range(len(mesh)):
#        rr[(mm == mesh[i])*(tt == time[i])] = r[i]


    ax.fill_between([0.8, 1.8], -1.5, 0.5, color='lightgrey')
    ax.pcolormesh(mm, ff, np.log10(rr), vmin=vmin, vmax=vmax, shading='nearest', rasterized=True)

    for i in range(len(M)):
        ax.text(M[i], FeH[i], '{:.2f}'.format(np.log10(r[i])), ha='center', va='center', size=9, color=sm_g.to_rgba(np.log10(r[i])))

    
    ax.set_xticks((0.9, 1.1, 1.3, 1.5, 1.7))
    ax.set_yticks((-1.2, -0.8, -0.4, 0.0, 0.4))
    ax.set_xlim(0.8, 1.8)
    ax.set_ylim(-1.5, 0.5)
    ax.invert_xaxis()
ax3.set_xlabel('Mass')
ax4.set_xlabel('Mass')
ax1.set_ylabel('[Fe/H]')
ax3.set_ylabel('[Fe/H]')


cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap('viridis'), norm=norm, orientation='horizontal', ticklocation='top')
cbar.set_label(r'$\log_{10}(r)$')

fig.savefig('mesa_r_spread.png', dpi=600, bbox_inches='tight')
fig.savefig('mesa_r_spread.pdf', dpi=600, bbox_inches='tight')


