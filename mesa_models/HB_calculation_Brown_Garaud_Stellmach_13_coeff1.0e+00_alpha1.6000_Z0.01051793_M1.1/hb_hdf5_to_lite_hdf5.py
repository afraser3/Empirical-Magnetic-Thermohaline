import pathlib
import numpy as np
import h5py

dense_hdf = 'hdf5_LOGS.h5'
lite_hdf  = 'hdf5_lite.h5'

here = pathlib.Path('./').resolve()
case = str(here).split('/')[-1].split('work_')[-1]

good_profiles = ['mass', 'kt', 'kmu', 'gradT_sub_grada', 'gradr_sub_grada', 'gradL_composition_term', \
                 'eps_nuc', 'mixing_type', 'radius', 'tau', 'zone',  'log_D_mix', \
                 'brunt_N2', 'brunt_N2_structure_term', 'logT', 'logRho', 'cp', 'opacity']

with h5py.File(dense_hdf, 'r') as in_f:
    with h5py.File(lite_hdf, 'w') as out_f:
        in_f.copy('history', out_f)

        #radial resolution changes over time so each profile is a different group
        profile_group = out_f.create_group('profiles')
        for i in [int(k) for k in in_f['profiles'].keys()]:
            print('checking profile {}'.format(i))
            eps_nuc = in_f['profiles/{}/bulk/eps_nuc'.format(i)][()]
            mixing_type = in_f['profiles/{}/bulk/mixing_type'.format(i)][()]
            mass = in_f['profiles/{}/bulk/mass'.format(i)][()]
            mass_coord = mass/mass.max()

            max_burn = np.argmax(eps_nuc)
            good_mass = (mass_coord < mass_coord[max_burn]*1.1)
            if max_burn != len(eps_nuc) - 1:
                good_mass *= (mass_coord > mass_coord[max_burn+1])
            else:
                good_mass *= (mass_coord > mass_coord[max_burn])
            thermohaline = mixing_type==4
            good_bool = good_mass*thermohaline
            if np.sum(good_bool):
                print('adding {} to lite file'.format(i))
                this_group = profile_group.create_group('{}'.format(i))
                in_f['profiles/{}'.format(i)].copy('headers', this_group)
                bulk = this_group.create_group('bulk')
                for k in good_profiles:
                    bulk[k] = in_f['profiles/{}/bulk/{}'.format(i,k)][()][good_mass]
