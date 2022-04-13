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


fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(6.5, 4))

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
    axs[k][j].scatter(np.log10(r), mixing_shet, marker='s', color=Dark2_3.mpl_colors[0], label='Shetrone', **kwargs)
    axs[k][j].scatter(np.log10(r), mixing_apk, marker='o', color=Dark2_3.mpl_colors[2], label='APOGEE', **kwargs)
    axs[k][j].text(0.98, 0.98, tag, ha='right', va='top', transform=axs[k][j].transAxes, fontsize=9)
    axs[k][j].set_xlim(-4.25, -2.75)
    axs[k][j].set_ylim(-0.05, 0.7)
    axs[k][j].set_xticks((-4.2, -3.8, -3.4, -3.0))

    if i == 0:
        axs[k][j].legend(loc='center right', fontsize=8)

axs[1][1].set_xlabel(r'$\log_{10} r$')
axs[1][0].set_xlabel(r'$\log_{10} r$')
axs[0][0].set_ylabel(r'$\Delta$ [C/N]')
axs[1][0].set_ylabel(r'$\Delta$ [C/N]')



fig.savefig('mixing_vs_r.png', dpi=300, bbox_inches='tight')
fig.savefig('mixing_vs_r.pdf', dpi=300, bbox_inches='tight')
