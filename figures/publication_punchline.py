import numpy as np
from scipy import optimize as opt
import matplotlib
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Dark2_4
plt.style.use('./apj.mplstyle')

from Nu_R0_model_comparison import plot_parameterizations

data = np.genfromtxt("./mixing_vs_r/scratchtableforEvan_figD.txt", skip_header=4)

M = data[:,0]
feh = data[:,1]
mixing_apk = data[:,2]
mixing_shet = data[:,3]
r = r_brown = data[:,4]

norm = matplotlib.colors.Normalize(vmin=feh.min(), vmax=feh.max())
sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
sm.set_array([])

fig = plt.figure(figsize=(6.5, 2))
ax1 = fig.add_axes((0.00, 0.00, 0.45, 1.00))
ax2 = fig.add_axes((0.55, 0.00, 0.45, 0.92))
cax = fig.add_axes((0.55, 0.95, 0.45, 0.04)) 

plot_parameterizations(ax1)


#fig 4 reproduction
kwargs = {'edgecolor' : 'k'}
for m in range(len(M)):
    ax2.scatter(r[m], mixing_shet[m], marker='s', color=sm.to_rgba(feh[m]), label='Shetrone', **kwargs)
    ax2.scatter(r[m], mixing_apk[m], marker='o', color=sm.to_rgba(feh[m]), label='APOKASC', **kwargs)
    if m == 0:
        ax2.legend(loc='center right', fontsize=8)
ax2.text(0.98, 0.98, 'Brown', ha='right', va='top', transform=ax2.transAxes, fontsize=9)
ax2.set_xlim(1e-5, 1)
ax2.set_ylim(-0.05, 0.7)
ax2.set_xscale('log')

cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap('plasma'), norm=norm, orientation='horizontal', ticklocation='top')
cbar.set_label(r'[Fe/H]')

ax1.set_xlabel(r'$r$')
ax2.set_xlabel(r'$r$')
ax2.set_ylabel(r'$\Delta$ [C/N]')


fig.savefig('punchline.png', dpi=300, bbox_inches='tight')
fig.savefig('punchline.pdf', dpi=300, bbox_inches='tight')




