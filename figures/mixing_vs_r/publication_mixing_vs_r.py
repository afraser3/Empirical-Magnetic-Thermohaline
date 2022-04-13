import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Dark2_3
plt.style.use('../apj.mplstyle')

data = np.genfromtxt("scratchtableforEvan_figD.txt", skip_header=4)

M = data[:,0]
feh = data[:,1]
mixing_apk = data[:,2]
mixing_shet = data[:,3]
r_brown = data[:,4]
r_kip_low = data[:,5]
r_kip_med = data[:,6]
r_kip_hi  = data[:,7]

norm = matplotlib.colors.Normalize(vmin=feh.min(), vmax=feh.max())
sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
sm.set_array([])




fig = plt.figure(figsize=(6.5, 4))
ax1 = fig.add_axes((0.00, 0.53, 0.45, 0.4))
ax2 = fig.add_axes((0.55, 0.53, 0.45, 0.4))
ax3 = fig.add_axes((0.00, 0.00, 0.45, 0.4))
ax4 = fig.add_axes((0.55, 0.00, 0.45, 0.4))
axs = [[ax1, ax2], [ax3, ax4]]

cax = fig.add_axes((0, 0.95, 1, 0.04)) 

for i in range(4):
    j = i % 2
    k = i // 2
    if i == 0:
        r = r_brown
        tag = 'Brown'
    elif i == 1:
        r = r_kip_low
        tag = r'Kipp, $\alpha_{\rm{th}} = 0.1$'
    elif i == 2:
        r = r_kip_med
        tag = r'Kipp, $\alpha_{\rm{th}} = 2$'
    elif i == 3:
        r = r_kip_hi
        tag = r'Kipp, $\alpha_{\rm{th}} = 700$'
    kwargs = {'edgecolor' : 'k'}
    for m in range(len(M)):
        #        axs[k][j].scatter(np.log10(r[m]), mixing_shet[m], marker='s', color=sm.to_rgba(feh[m]), label='Shetrone', **kwargs)
        #        axs[k][j].scatter(np.log10(r[m]), mixing_apk[m], marker='o', color=sm.to_rgba(feh[m]), label='APOGEE', **kwargs)
        axs[k][j].scatter(r[m], mixing_shet[m], marker='s', color=sm.to_rgba(feh[m]), label='Shetrone', **kwargs)
        axs[k][j].scatter(r[m], mixing_apk[m], marker='o', color=sm.to_rgba(feh[m]), label='APOGEE', **kwargs)
        if i == 0 and m == 0:
            axs[k][j].legend(loc='center right', fontsize=8)
    axs[k][j].text(0.98, 0.98, tag, ha='right', va='top', transform=axs[k][j].transAxes, fontsize=9)
    axs[k][j].set_xlim(10**(-4.25), 10**(-2.75))
    axs[k][j].set_ylim(-0.05, 0.7)
#    axs[k][j].set_xticks((1e-4, 3e-4, 1e-3))
    axs[k][j].set_xscale('log')
#    axs[k][j].set_xticks((-4.2, -3.8, -3.4, -3.0))


cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap('plasma'), norm=norm, orientation='horizontal', ticklocation='top')
cbar.set_label(r'[Fe/H]')

axs[1][1].set_xlabel(r'$r$')
axs[1][0].set_xlabel(r'$r$')
axs[0][0].set_ylabel(r'$\Delta$ [C/N]')
axs[1][0].set_ylabel(r'$\Delta$ [C/N]')



fig.savefig('mixing_vs_r.png', dpi=300, bbox_inches='tight')
fig.savefig('mixing_vs_r.pdf', dpi=300, bbox_inches='tight')
