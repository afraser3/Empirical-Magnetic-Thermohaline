import pathlib
import mesa_reader as mr
import h5py

here = pathlib.Path('./').resolve()
case = str(here).split('/')[-1].split('work_')[-1]
logdirkw = dict()
logdirkw['profile_prefix'] = '{}_no-rot-minDmix.profile'.format(case)
logdirkw['index_file'] = '{}s.index'.format(logdirkw['profile_prefix'])
logdirkw['history_file'] = 'history_{}_no-rot-minDmix.data'.format(case.split('thermohaline_')[-1])

l = mr.MesaLogDir(**logdirkw)
for k in dir(l):
    print(k)

with h5py.File('hdf5_LOGS.h5', 'w') as f:
    history_group = f.create_group('history')
    history_headers = history_group.create_group('headers')
    history_bulk    = history_group.create_group('bulk')

    for k in l.history_data.header_names:
        history_headers[k] = l.history_data.header_data[k]
    for k in l.history_data.bulk_names:
        history_bulk[k] = l.history_data.bulk_data[k]

    #radial resolution changes over time so each profile is a different group
    profile_group = f.create_group('profiles')
    for i in l.profile_numbers:
        print('stitching profile {}'.format(i))
        this_group = profile_group.create_group('{}'.format(i))
        headers = this_group.create_group('headers')
        bulk = this_group.create_group('bulk')
        data = l.profile_data(profile_number=i)

        for k in data.header_names:
            headers[k] = data.header_data[k]
        for k in data.bulk_names:
            bulk[k] = data.data(k)
