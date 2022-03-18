#!/usr/bin/env python3
####################################################
#
# Author: M Joyce
#
####################################################
import numpy as np
import glob
import sys
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import argparse

import h5py


tag = 'grid2'

colors=['red', 'orange','lightgreen',  'green', 'blue',\
        'navy','indigo', 'purple', 'violet',\
        'magenta', 'pink', 'maroon', 'brown', 'black','gray',\
        'red', 'orange','lightgreen',  'green', 'blue',\
        'navy','indigo', 'purple', 'violet',\
        'magenta', 'pink', 'maroon', 'brown', 'black','gray',\
        'red', 'orange','lightgreen',  'green', 'blue',\
        'navy','indigo', 'purple', 'violet',\
        'magenta', 'pink', 'maroon', 'brown', 'black','gray']

fig, ax = plt.subplots(figsize=(18,16))

with h5py.File('hdf5_LOGS.h5', 'r') as f:
    history = f['history/bulk']

    mesa_histories = dict()
    for key in ['model_number', 'star_mass', 'star_age', 'log_LH', 'log_Teff', 'log_L', 'log_R', 'log_g']:
        mesa_histories[key] = history[key][()]
    locals().update(mesa_histories)

    Teff = 10.0**log_Teff
    star_age = star_age/1e9
    radius = 10.0**log_R


    ax.plot(Teff, log_L, '-', lw=4, ms=1, color=colors[0])
    


ax.set_xlabel(r"Teff (K)",  fontsize=36)
ax.set_ylabel(r"Log $L/L_{\odot}$", fontsize=36)

ax.tick_params(axis='both', which='major', labelsize=24)
ax.tick_params(axis='both', which='minor', labelsize=20)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.tick_params(which='both', width=4)
ax.tick_params(which='major', length=12)
ax.tick_params(which='minor', length=8, color='black')

#plt.title(r'Adjustments to $\alpha_{MLT}$ only (No $B$, no spots)', fontsize = 30)

#plt.xlim(3280,4000)
#plt.ylim(-0.2,0.41)
plt.gca().invert_xaxis()
plt.legend(loc=6, fontsize=16)

plt.savefig('HR_diagram.png')
plt.close()
