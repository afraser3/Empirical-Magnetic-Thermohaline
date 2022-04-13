import pathlib
import glob

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


csvfiles = glob.glob('*.csv')
for fname in csvfiles:
    vmin = -4
    vmax = -3
    name = fname.split('.csv')[0]

    data = np.genfromtxt(fname, delimiter=',', skip_header=1)
    M = data[:,0]
    FeH = data[:,1]
    r = data[:,3]
#    mesh = data[:,6]
#    time = data[:,7]

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


    fig = plt.figure()
    ax = fig.add_axes((0.1, 0.1, 0.8, 0.7))
    cax = fig.add_axes((0.1, 0.94, 0.8, 0.03))

#    ax.pcolormesh(tt, mm, np.log10(rr), vmin=vmin, vmax=vmax, shading='nearest')
    ax.pcolormesh(mm, ff, np.log10(rr), vmin=vmin, vmax=vmax, shading='nearest')

    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap='viridis', norm=norm)
    sm.set_array([])

    for i in range(len(M)):
        ax.text(M[i], FeH[i], '{:.3f}'.format(np.log10(r[i])), ha='center', va='center')


    ax.set_xlabel('Mass')
    ax.set_ylabel('Fe/H')
    ax.invert_xaxis()
#    ax.set_xlabel('mesh delta coeff')
#    ax.set_ylabel('time delta coeff')
#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    ax.invert_yaxis()



    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap('viridis'), norm=norm, orientation='horizontal')
    cbar.set_label('r')

    fig.savefig('convergence_{}.png'.format(name), dpi=300, bbox_inches='tight')
