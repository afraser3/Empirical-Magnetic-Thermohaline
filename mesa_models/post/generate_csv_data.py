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



z_map = np.around(np.genfromtxt('../Z_feH_table.dat', skip_header=1, usecols=(0, 3)), decimals=8)
FeHvals   = z_map[:,0].ravel()
Zvals = z_map[:,1].ravel()

here = pathlib.Path('./')
def read_dirs(model_name, coeff):
    data = []
    for path in natural_sort(here.glob('work*/')):
        for strval in str(path).split('_'):
            if 'M' in strval:
                M = float(strval.split('M')[-1])
            if 'Z' in strval:
                Z = float(strval.split('Z')[-1])
                FeH = FeHvals[Zvals == Z][0]
        try:
            if model_name not in str(path):
                continue
            if float(str(path).split('coeff')[-1].split('_')[0]) != coeff:
                continue
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
#                N = int(np.sum(good_range))
#                rval = np.median(r[good_range][N//2:])
#                R0val = np.median(R0[good_range][N//2:])
#                tauval = np.median(tau[good_range][N//2:])
                rval = np.median(r[good_range])
                R0val = np.median(R0[good_range])
                tauval = np.median(tau[good_range])
                plt.xlabel('model number')
                plt.ylabel('r')

                plt.savefig(str(path.joinpath('r_vs_num.png')), dpi=400)
                plt.clf()
                data.append((M, FeH, Z, rval, R0val, tauval))
                print(data[-1])
        except:
            print('ERROR: r_vs_time.h5 not found in {}'.format(path))
            data.append((M, FeH, Z, np.nan, np.nan, np.nan))
    if np.array(data).size > 0:
#        data = sorted(data, key=lambda d: (d[0], d[1]))
        for i in range(len(data)):
            print(data[i])
        header = "{:>18s}".format("M") + (5*"{:>21s}").format("Fe/H", "Z", "r", "R0", "tau")
        np.savetxt('mesa_{}_coeff{:.1e}_values.csv'.format(model_name, coeff), data, '%20.10f', delimiter=',', header=header)


read_dirs('Brown', 1)
read_dirs('Kipp', 0.1)
read_dirs('Kipp', 2)
read_dirs('Kipp', 700)
