import pathlib
import re

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import mode

def natural_sort(iterable, reverse=False):
    """
    Sort alphanumeric strings naturally, i.e. with "1" before "10".
    Based on http://stackoverflow.com/a/4836734.

    Taken from dedalus/tools/general.py

    """

    convert = lambda sub: int(sub) if sub.isdigit() else sub.lower()
    key = lambda item: [convert(sub) for sub in re.split('([0-9]+)', str(item))]

    return sorted(iterable, key=key, reverse=reverse)



z_map = np.around(np.genfromtxt('../../Z_feH_table.dat', skip_header=1, usecols=(0, 3)), decimals=8)
FeHvals   = z_map[:,0].ravel()
Zvals = z_map[:,1].ravel()

here = pathlib.Path('./')
def read_dirs(massval, zval):
    data = []
    for path in natural_sort(here.glob('work*/')):
        for strval in str(path).split('_'):
            if 'M' in strval:
                M = float(strval.split('M')[-1])
            if 'Z' in strval:
                Z = float(strval.split('Z')[-1])
                FeH = FeHvals[Zvals == Z][0]
            if 'mesh' in strval:
                mesh = float(strval.split('mesh')[-1])
            if 'time' in strval:
                time = float(strval.split('time')[-1])
        if M != massval or Z != zval:
            continue
        try:
            with h5py.File(path.joinpath('r_vs_time.h5'), 'r') as f:
#                print('reading {}'.format(str(path)))

                model = f['models'][()]
                shortmodel, unique = np.unique(model, return_index=True)
                model = model[unique]

                good_range = f['good_point'][()][unique]
                if np.sum(good_range) == 0: continue
                r = f['r'][()][unique]
                R0 = f['R0'][()][unique]
                tau = f['tau'][()][unique]
                good_r = r[good_range]
                good_model = model[good_range]
                plt.axhline(np.median(good_r))
                plt.plot(model, r)
                plt.plot(good_model, good_r)
                plt.ylim(1e-4, 1)
                plt.yscale('log')
                rval = np.median(r[good_range])
                R0val = np.median(R0[good_range])
                tauval = np.median(tau[good_range])
                plt.xlabel('model number')
                plt.ylabel('r')

                plt.savefig(str(path.joinpath('r_vs_num.png')), dpi=400)
                plt.clf()
                data.append((M, FeH, Z, rval, R0val, tauval, mesh, time))
        except:
            print('ERROR: r_vs_time.h5 not found in {}'.format(path))
            data.append((M, FeH, Z, np.nan, np.nan, np.nan, np.nan, np.nan))
    if np.array(data).size > 0:
#        for i in range(len(data)):
#            print(data[i])
        header = "{:>18s}".format("M") + (7*"{:>21s}").format("Fe/H", "Z", "r", "R0", "tau", "mesh_coeff", "time_coeff")
        M, FeH = data[0][:2]
        np.savetxt('mesa_FeH{:.1f}_Mv{:.1f}_values.csv'.format(FeH, M), data, '%20.8f', delimiter=',', header=header)
        print('saved ' + 'mesa_FeH{:.1f}_Mv{:.1f}_values.csv'.format(FeH, M))


masses = [0.9, 1.3, 1.7]
Zs = [0.03814436, 0.00671765, 0.00108379]

for Zv in Zs:
    for Mv in masses:
        read_dirs(Mv, Zv)
