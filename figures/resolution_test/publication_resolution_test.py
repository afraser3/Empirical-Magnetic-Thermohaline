import pathlib
import glob

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

plt.style.use('../apj.mplstyle')

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

axs = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]

ref_file = '../mesa_spread/mesa_Brown_coeff1.0e+00_values.csv'
ref_data = np.genfromtxt(ref_file, delimiter=',', skip_header=1)
ref_Z = ref_data[:,2]
ref_M = ref_data[:,0]
ref_r = ref_data[:,3]

csvfiles = glob.glob('*.csv')
err_cmap = 'Reds'
vmin = -3
vmax = 2
norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
sm = plt.cm.ScalarMappable(cmap=err_cmap, norm=norm)
sm.set_array([])
#
#norm_g = matplotlib.colors.Normalize(vmin=vmin + 0.5, vmax=vmax)
#sm_g = plt.cm.ScalarMappable(cmap='Greys', norm=norm)
#sm_g.set_array([])
#
for fname in csvfiles:
    Z = float(fname.split('Z')[-1].split('_')[0])
    M = float(fname.split('M')[-1].split('.')[0])
    data = np.genfromtxt(fname, delimiter=',', skip_header=1)

    grid_FeH = ref_data[:,1]
    grid_Z = ref_data[:,2]
    grid_M = ref_data[:,0]
    mesh = data[:,6]
    time = data[:,7]
    r = data[:,3]

    

    #hack
    if Z == 0.028:
        Z = 0.02517204
        FeH = 0.2
    if M == 1:
        M = 1.1
    reference_r = r[(mesh == 0.5)*(time == 0.1)]
    print(Z, M)
    err = np.abs(1 - r/reference_r)

    #Fix this once real data exist
#    FeH = grid_FeH[grid_Z == Z][0]
#    err = np.abs(1 - r/ref_r[(ref_Z == Z)*(ref_M == M)])


    mm, tt = np.meshgrid(np.unique(mesh), np.unique(time))
    cc = np.zeros_like(tt)
    cc[:] = np.nan
    for i in range(len(time)):
        cc[(mm == mesh[i])*(tt == time[i])] = err[i]
        if cc[(mm == mesh[i])*(tt == time[i])] == 0:
            cc[(mm == mesh[i])*(tt == time[i])] = 1e-8
    print(cc)
    for ax in axs:
        ax.fill_between([0, 2], 0, 2, color='lightgrey')
        ax.pcolormesh(tt, mm, np.log10(cc), vmin=vmin, vmax=vmax, shading='nearest', rasterized=True, cmap='Reds')

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlim(1e-1, 1.05e0)
        ax.set_ylim(1e-1, 1.05e0)

        bbox_props = dict(ec="k", alpha=1, fc="w", linewidth=0.8, pad=2)
        ax.text(0.017, 1.055, 'M = {:.1f}, [Fe/H] = {:.1f}'.format(M, FeH), transform=ax.transAxes, ha='left', va='center', size=9, bbox=bbox_props)
        

ax1.set_ylabel('mesh coeff')
ax4.set_ylabel('mesh coeff')
ax7.set_ylabel('mesh coeff')

ax7.set_xlabel('time coeff')
ax8.set_xlabel('time coeff')
ax9.set_xlabel('time coeff')


cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap(err_cmap), norm=norm, orientation='horizontal', ticklocation='top')
cbar.set_label(r'fractional error')

fig.savefig('resolution_test.png', dpi=600, bbox_inches='tight')
fig.savefig('resolution_test.pdf', dpi=600, bbox_inches='tight')


