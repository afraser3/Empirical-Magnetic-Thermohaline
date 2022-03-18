####################################################
#
# Author: M Joyce & E Anders
# Use on command line as python3 get_fluid_parameters.py
# Output is stored in fluid_parameter_plots/ folder
#
####################################################
import pathlib
import os
import numpy as np
import glob
import sys
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import argparse
import matplotlib.patches as patches
import h5py


#This is a bad idea but the divide by zero errors drive me insane
import warnings
warnings.filterwarnings("ignore")

#---------------------------------
savefig = True
colors=['red', 'orange','lightgreen',  'green', 'blue',\
        'navy','indigo', 'purple', 'violet',\
        'magenta', 'pink', 'maroon', 'brown', 'black','gray']

def compute_R0(gradT_sub_grada, gradr_sub_grada, gradL_composition_term, *args, **kwargs):
    use_gradT = kwargs.get('use_gradT', True)
    delta =  1 ## note: not -1 because of order of differences
    phi = 1
    
    if use_gradT:
        R0 = (delta/phi)*(gradT_sub_grada)/(gradL_composition_term)
    else:
        R0 = (delta/phi)*(gradr_sub_grada)/(gradL_composition_term)
        print('using gradr_sub_grada')
    return R0

def compute_tau(kmu, kt):
    tau = kmu/kt
    return tau

#################################################
#
# 501:   !mixing_type ! mixing types are defined in mesa/const/public/const_def  
#
# to identify mixing type:
#
#       ! mixing types
#       ! NOTE: some packages may depend on the order
#       integer, parameter :: no_mixing = 0
#       integer, parameter :: convective_mixing = 1
#       integer, parameter :: overshoot_mixing = 2
#       integer, parameter :: semiconvective_mixing = 3
#       integer, parameter :: thermohaline_mixing = 4
#       integer, parameter :: rotation_mixing = 5
#       integer, parameter :: rayleigh_taylor_mixing = 6
#       integer, parameter :: minimum_mixing = 7
#       integer, parameter :: anonymous_mixing = 8  ! AKA "WTF_mixing"
#       integer, parameter :: leftover_convective_mixing = 9  ! for regions with non-zero conv_vel that are not unstable to convection
#                                                             ! used for time dependent convection
#     
#       integer, parameter :: number_of_mixing_types = leftover_convective_mixing+1
#
#################################################

import re
def natural_sort(iterable, reverse=False):
    """
    Sort alphanumeric strings naturally, i.e. with "1" before "10".
    Based on http://stackoverflow.com/a/4836734.

    Natural sort function from dedalus/tools/general.py
    See https://github.com/DedalusProject/dedalus/blob/master/dedalus/tools/general.py
    """

    convert = lambda sub: int(sub) if sub.isdigit() else sub.lower()
    key = lambda item: [convert(sub) for sub in re.split('([0-9]+)', str(item))]

    return sorted(iterable, key=key, reverse=reverse)

output_folder = 'fluid_parameter_plots/'
if not os.path.exists(output_folder):
    os.mkdir(output_folder)

fig, ax = plt.subplots(figsize=(14,10))
with h5py.File('hdf5_LOGS.h5', 'r') as f:
    here = pathlib.Path('./').resolve()
    here_str = str(here).split('/')[-1]
    mass_str = here_str.split('_M')[1].split('_')[0]
    alpha = here_str.split('_alpha')[1].split('_Z')[0]
    Z = here_str.split('_Z')[1].split('_M')[0]
    for i in range(len(f['profiles'].keys())):
        lognum = i+1
        bulk_group = f['profiles'][str(lognum)]['bulk']
        model_number = f['history']['bulk']['model_number'][()][i]
        print('reading profile {} / model number {}'.format(lognum, model_number))

        mesa_profiles = dict()
        for key in ['mass', 'kt', 'kmu', 'gradT_sub_grada', 'gradr_sub_grada', 'gradL_composition_term',\
                    'eps_nuc', 'mixing_type', 'radius', 'tau', 'pressure_scale_height', 'log_D_mix', 'logRho',\
                    'Prandtl', 'logT' ]:
            mesa_profiles[key] = bulk_group[key][()]
        locals().update(mesa_profiles)

        th=np.where((mixing_type==4)) ## worry about this mass edge?    
        try:
            kt = 10**kt
            kmu = 10**kmu
            e_nuc = eps_nuc
            # vis = mesa.vis
            # tau = kmu/kt 


            R0 = compute_R0(gradT_sub_grada, gradr_sub_grada, gradL_composition_term, use_gradT = True)
            tau = compute_tau(kmu, kt)
            mass_coord = mass/mass.max()


            plt.plot(mass_coord, logRho,\
             linewidth=3, color='green', label=r'Log $\rho$')
            plt.plot(mass_coord, np.log10(pressure_scale_height),\
             linewidth=3, color='orange', label=r'log $H_p$')
            plt.plot(mass_coord, Prandtl,\
             linewidth=3, linestyle = '-',  color='yellow', label=r'Pr')
            plt.plot(mass_coord, np.log10(R0),\
             linewidth=3, linestyle = '-',  color='palevioletred', label=r'log $R_0$')
            plt.plot(mass_coord, np.log10(tau),\
             linewidth=3, linestyle = ':',  color='cornflowerblue', label=r'log $\tau =$ log $k_{\mu}/k_{T}$')
            plt.plot(mass_coord, np.log10(e_nuc),\
             linewidth=3, linestyle = '-',  color='black', label=r'log $\varepsilon_{nuc}$')
            plt.plot(mass_coord, log_D_mix/3.0,\
             linewidth=3, linestyle = '--',  color='blue', label=r'log Dmix/3 ')
            plt.plot(mass_coord, logT/2.0,\
             linewidth=3, linestyle = '--',  color='orange', label=r'(logT)/2')

            plt.ylabel('Various', fontsize = 24)


            ## Imported from Matteo's scripts on GitHub
            #th=np.where((mixing_type==4) & (mass_coord > 0.2) & (mass_coord < 0.28)) ## worry about this mass edge?
            #R0_region_of_interest = np.where( (R0 > 1.0) & ( R0 < 1.0/tau) )[0]
            R0_region_of_interest = th
            min_f0= mass_coord[R0_region_of_interest].min() #np.log10(3e6)#x
            min_f1= -7.8 #0.5e9 #y
            max_f0= mass_coord[R0_region_of_interest].max()
            max_f1= 7.8 #2e9 #log_max_T.max()

            # R0_region_of_interest = np.where( (R0 < 1.0) & ( R0 > 1.0/tau) )[0]
            # min_f0= mass_coord[R0_region_of_interest].max() #np.log10(3e6)#x
            # min_f1= -7.8 #0.5e9 #y
            # max_f0= mass_coord[R0_region_of_interest].min()
            # max_f1= 7.8 #2e9 #log_max_T.max()


            width = max_f0 - min_f0
            height = max_f1 - min_f1
            ax.add_patch(
                patches.Rectangle(
                    xy=(min_f0, min_f1),  # point of origin.
                    width=width,
                    height=height,
                    linewidth=1,
                    color='grey',
                    alpha=0.3,
                    fill=True,
                    label = r'Mixing_Type = thermohaline'
                    #label=r'$1 < R_0 < 1/\tau$'
                )
            )

            plt.legend(loc=1, fontsize = 16)

            ax.set_xlabel(r"Mass/Mstar",  fontsize=36)
            ax.set_ylabel(r"Various", fontsize=36)

            ax.tick_params(axis='both', which='major', labelsize=24)
            ax.tick_params(axis='both', which='minor', labelsize=20)

            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(which='both', width=4)
            ax.tick_params(which='major', length=12)
            ax.tick_params(which='minor', length=8, color='black')

            textstr=r'm='+str(mass_str)+'\nZ='+str(Z) +'\nprofile='+str(lognum)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            ax.text(0.305, 0, textstr, fontsize=16, bbox = props)


            if savefig:
                try:
                    plt.savefig('{}/fluid_quantities_m{}_Z{}_profile{:06d}.png'.format(output_folder, mass_str, Z, lognum))
                except:
                    plt.savefig('{}/fluid_quantities_m{}_Z{}.png'.format(output_folder, mass_str, Z))
            else:
                plt.show()
            plt.clf()


        except ValueError:
            plt.clf()
            print("")
            print("Thermohaline mixing not detected in model m=" + str(mass_str) + '; Z='+str(Z) + '; profile_num = {:06d}'.format(lognum) )
            print("")


