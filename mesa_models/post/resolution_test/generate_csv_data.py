import pathlib

import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import mode


z_map = np.around(np.genfromtxt('../../Z_feH_table.dat', skip_header=1, usecols=(0, 3)), decimals=4)
FeHvals   = z_map[:,0].ravel()
Zvals = z_map[:,1].ravel()

print(Zvals)

here = pathlib.Path('./')
def read_dirs(model_name, coeff):
    data = []
    for path in natural_sort(here.glob('work*/')):
        try:
            if model_name not in str(path):
                continue
            if float(str(path).split('coeff')[-1].split('_')[0]) != coeff:
                continue
            with h5py.File(path.joinpath('r_vs_time.h5'), 'r') as f:
#                print('reading {}'.format(str(path)))
                for strval in str(path).split('_'):
                    if 'M' in strval:
                        M = float(strval.split('M')[-1])
                    if 'Z' in strval:
                        Z = float(strval.split('Z')[-1])
                        FeH = FeHvals[Zvals == Z][0]
                    if 'mesh' in strval:
                        mesh = float(strval.split('mesh')[-1])
                    if 'time' in strval:
                        time = float(strval.split('time')[-1]))


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
                print(data[-1])
        except:
            print('ERROR: r_vs_time.h5 not found in {}'.format(path))
    if np.array(data).size > 0:
        for i in range(len(data)):
            print(data[i])
        header = "{:>18s}".format("M") + (7*"{:>21s}").format("Fe/H", "Z", "r", "R0", "tau", "mesh", "time")
        np.savetxt('mesa_{}_coeff{:.1e}_values.csv'.format(model_name, coeff), data, '%20.8f', delimiter=',', header=header)







read_dirs('Brown', 1)
#read_dirs('Kipp', 0.1)
#read_dirs('Kipp', 2)
#read_dirs('Kipp', 80)
#read_dirs('Kipp', 700)
