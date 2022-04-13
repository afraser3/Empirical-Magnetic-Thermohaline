"""
This python script generates MESA run directories based on a template directory 
for a specified range of metallicities and initial masses.
"""
import shutil
from pathlib import Path
import numpy as np

Z_values = [0.00108379,  0.00671765,  0.03814436] #fe/h = -1.2, -0.4, 0.4
M_values = [0.9, 1.3, 1.7 ]
mesh_values = np.logspace(-1, 0, 5)
time_values = np.logspace(-1, 0, 5)
th_model = ["'Brown_Garaud_Stellmach_13'",]
th_coeff = [1]
alpha_MLT = 1.6
a_string  = '{:.2f}'.format(alpha_MLT)


#Z - to - Y from Jamie: https://ui.adsabs.harvard.edu/abs/2022ApJ...927...31T/abstract
z_solar = 0.0179492
y_solar = 0.2725693
delta_y_div_delta_z = 1.3426 
big_bang_helium = 0.2485 
z_to_y = lambda zval: big_bang_helium + delta_y_div_delta_z * zval

template_dir = Path('../template_work_thermohaline/')

for thermo_model, thermo_coeff in zip(th_model, th_coeff):
    for mesh_delta in mesh_values:
        for time_delta in time_values:
            for Zv in Z_values:
                Yv = z_to_y(Zv)
                for Mv in M_values:
                    zv_string = '{:.8f}'.format(Zv)
                    yv_string = '{:.8f}'.format(Yv)
                    mv_string = '{:.2f}'.format(Mv)
                    mesh_string = '{:.6f}'.format(mesh_delta)
                    time_string = '{:.6f}'.format(time_delta)
                    run_tag = 'mesh{}_time{}_alpha{}_Z{}_M{}'.format(mesh_string, time_string, a_string, zv_string, mv_string)
                    run_dir = Path('./work_thermohaline_{}'.format(run_tag))

                    if not run_dir.exists():
                        print('creating directory {}'.format(run_dir))
                        run_dir.mkdir()
                    else:
                        print('{} already exists'.format(run_dir))

                    for filename in ['clean', 'mk', 're', 'rn', 'inlist', 'inlist_pgstar_thermohaline', 'controls.defaults', 'history_columns_thermohaline.list', 'profile_columns_thermohaline.list', 'star_job.defaults', 'mesa_to_hdf5.py', 'get_fluid_parameters.py', 'find_r_vs_time.py']:
                        template_path = template_dir.joinpath(filename)
                        run_path = run_dir.joinpath(filename)
                        shutil.copy(str(template_path), str(run_path))

                    for directory in ['make', 'src']:
                        template_path = template_dir.joinpath(directory)
                        run_path = run_dir.joinpath(directory)
                        if not run_path.exists():
                            shutil.copytree(str(template_path), str(run_path))

                    #copy inlist_project_temp w/ appropriate changes
                    template_inlist = open(str(template_dir.joinpath('inlist_project_temp')), 'r')
                    run_inlist = open(str(run_dir.joinpath('inlist_project_temp')), 'w')
                    replacements = dict()
                    mod_run_tag = '{}_no-rot-minDmix'.format(run_tag)
                    replacements['Zbase'] = zv_string
                    replacements['initial_z'] = zv_string
                    replacements['initial_y'] = yv_string
                    replacements['initial_mass'] = mv_string
                    replacements['save_model_filename'] = "'thermohaline_{}.mod'".format(mod_run_tag)
                    replacements['filename_for_profile_when_terminate'] = "'thermohaline_{}.profile.data'".format(mod_run_tag)
                    replacements['star_history_name'] = "'history_{}.data'".format(mod_run_tag)
                    replacements['profile_data_prefix'] = "'thermohaline_{}.profile'".format(mod_run_tag)
                    replacements['profiles_index_name'] = "'thermohaline_{}.profiles.index'".format(mod_run_tag)
#                    replacements['x_ctrl(12)'] = '1.25'
#                    replacements['x_ctrl(13)'] = '3.00'
                    replacements['x_ctrl(15)'] = mesh_delta
                    replacements['x_ctrl(16)'] = time_delta

                    replacements['thermohaline_option'] = thermo_model
                    replacements['thermohaline_coeff'] = '{:.10f}'.format(thermo_coeff)


#                    if Zv == 0.0280 and Mv == 1.7:
#                        replacements['x_ctrl(20)'] = '2.005d9'
#                        replacements['x_ctrl(21)'] = '2.030d9'
                    for line in template_inlist.readlines():
                        for key, string in replacements.items():
                            if key in line and 'initial_zfracs' not in line:
                                line = '{} = {}\n'.format(line.split('=', 1)[0], string)
                        run_inlist.write(line)
                    template_inlist.close()
                    run_inlist.close()

                    #copy evan's pleiades job submission file (may need to change for other people)
                    template_job = open(str(template_dir.joinpath('evan_pleiades_job_submit')), 'r')
                    run_job = open(str(run_dir.joinpath('evan_pleiades_job_submit')), 'w')
                    for line in template_job.readlines():
                        if "#PBS -N mesa_thermohaline" in line:
                            line = "#PBS -N mesa_thermohaline_{}\n".format(run_tag)
                        run_job.write(line)
                    template_job.close()
                    run_job.close()




