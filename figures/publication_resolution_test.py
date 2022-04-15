import pathlib
import glob
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def natural_sort(iterable, reverse=False):
    """
    Sort alphanumeric strings naturally, i.e. with "1" before "10".
    Based on http://stackoverflow.com/a/4836734.

    Taken from dedalus/tools/general.py
    """

    convert = lambda sub: int(sub) if sub.isdigit() else sub.lower()
    key = lambda item: [convert(sub) for sub in re.split('([0-9]+)', str(item))]

    return sorted(iterable, key=key, reverse=reverse)

plt.style.use('./apj.mplstyle')

fig = plt.figure(figsize=(6.5, 6.5))
ax1  = fig.add_axes((0.00, 0.64, 0.28, 0.25))
ax2  = fig.add_axes((0.36, 0.64, 0.28, 0.25))
ax3  = fig.add_axes((0.72, 0.64, 0.28, 0.25))
ax4  = fig.add_axes((0.00, 0.32, 0.28, 0.25))
ax5  = fig.add_axes((0.36, 0.32, 0.28, 0.25))
ax6  = fig.add_axes((0.72, 0.32, 0.28, 0.25))
ax7  = fig.add_axes((0.00, 0.00, 0.28, 0.25))
ax8  = fig.add_axes((0.36, 0.00, 0.28, 0.25))
ax9  = fig.add_axes((0.72, 0.00, 0.28, 0.25))
cax = fig.add_axes((0, 0.96, 1, 0.04))

axs = [ax3, ax2, ax1, ax6, ax5, ax4, ax9, ax8, ax7]

ref_file = './mesa_spread/mesa_Brown_coeff1.0e+00_values.csv'
ref_data = np.genfromtxt(ref_file, delimiter=',', skip_header=1)
ref_FeH = ref_data[:,1]
ref_Z = ref_data[:,2]
ref_M = ref_data[:,0]
ref_r = ref_data[:,3]

csvfiles = natural_sort(glob.glob('resolution_test/mesa*.csv'))
err_cmap = 'PiYG_r'
vmin = 0
vmax = 2
norm = matplotlib.colors.TwoSlopeNorm(np.log10(5), vmin=vmin-0.5, vmax=vmax)
#norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=err_cmap, norm=norm)
sm.set_array([])
#
#norm_g = matplotlib.colors.Normalize(vmin=vmin + 0.5, vmax=vmax)
#sm_g = plt.cm.ScalarMappable(cmap='Greys', norm=norm)
#sm_g.set_array([])
#
for i, fname in enumerate(csvfiles):
    print(fname)
    FeH = float(fname.split('FeH')[-1].split('_')[0])
    M = float(fname.split('Mv')[-1].split('_')[0])
    data = np.genfromtxt(fname, delimiter=',', skip_header=1)

    grid_FeH = ref_data[:,1]
    grid_Z = ref_data[:,2]
    grid_M = ref_data[:,0]
    mesh = data[:,6]
    time = data[:,7]
    r = data[:,3]

    #hack
    print(ref_FeH, ref_M, FeH, M)
    reference_r = ref_r[(ref_FeH == FeH)*(ref_M == M)]
#    r[(mesh == 0.5)*(time == 0.1)]
    err = 100*np.abs(1 - r/reference_r)

    #Fix this once real data exist
#    FeH = grid_FeH[grid_Z == Z][0]
#    err = np.abs(1 - r/ref_r[(ref_Z == Z)*(ref_M == M)])


    mm, tt = np.meshgrid(np.unique(mesh[np.isfinite(mesh)]), np.unique(time[np.isfinite(time)]))
    cc = np.zeros_like(tt)
    cc[:] = np.nan
    for j in range(len(time)):
        cc[(mm == mesh[j])*(tt == time[j])] = err[j]
        if cc[(mm == mesh[j])*(tt == time[j])] < 10**(vmin):
            cc[(mm == mesh[j])*(tt == time[j])] = 10**(vmin)
    print(cc)
    ax = axs[i]
    ax.fill_between([-1.1, 0.1], -1.1, 0.1, color='grey')
    ax.pcolormesh(np.log10(tt), np.log10(mm), np.log10(cc), norm=norm, shading='nearest', rasterized=True, cmap=err_cmap)

#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    ax.set_xlim(1e-1, 1.05e0)
#    ax.set_ylim(1e-1, 1.05e0)

    ax.set_xlim(-1.1, 0.1)
    ax.set_ylim(-1.1, 0.1)

    bbox_props = dict(ec="k", alpha=1, fc="w", linewidth=0.8, pad=2)
    ax.text(0.017, 1.055, 'M = {:.1f}, [Fe/H] = {:.1f}'.format(M, FeH), transform=ax.transAxes, ha='left', va='center', size=9, bbox=bbox_props)
    ax.scatter(-1, np.log10(0.5), c='k', marker='*')

    ax.set_xticks((-1, -0.5, 0))
    ax.set_yticks((-1, -0.5, 0))
    ax.set_xticklabels((r'$10^{-1}$', r'$10^{-0.5}$', r'$1$'))
    ax.set_yticklabels((r'$10^{-1}$', r'$10^{-0.5}$', r'$1$'))
        

ax1.set_ylabel('mesh coeff')
ax4.set_ylabel('mesh coeff')
ax7.set_ylabel('mesh coeff')

ax7.set_xlabel('time coeff')
ax8.set_xlabel('time coeff')
ax9.set_xlabel('time coeff')


cax.set_xlim(vmin, vmax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap(err_cmap), norm=norm, orientation='horizontal', ticklocation='top', boundaries=np.linspace(vmin, vmax, 1000))
cbar.set_label(r'percent difference')
cbar.set_ticks((0, np.log10(2.5), np.log10(5), 1, np.log10(25), np.log10(50), 2))
cbar.set_ticklabels((1, 2.5, 5, 10, 25, 50, 100))

fig.savefig('resolution_test.png', dpi=600, bbox_inches='tight')
fig.savefig('resolution_test.pdf', dpi=600, bbox_inches='tight')


