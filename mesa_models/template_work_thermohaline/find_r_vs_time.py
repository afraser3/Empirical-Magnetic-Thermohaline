import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
import palettable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

PLOT_MOVIE=True

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

if not os.path.exists('./figs_R0/'):
    os.mkdir('./figs_R0/')

good_models = []
good_ages = []
good_R0s = []
good_taus = []
good_rs = []
good_log_g = []

drs = []
dms = []
tau_mixs = []
good_point = []

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
axins = inset_axes(ax, width=1.3, height=0.9)
with h5py.File('hdf5_lite.h5', 'r') as f:
    log_gs = f['history']['bulk']['log_g'][()]
    models = f['history']['bulk']['model_number'][()]
    star_age = f['history']['bulk']['star_age'][()]
    star_mass = f['history']['bulk']['star_mass'][()]

    for lognum in sorted([int(k) for k in f['profiles'].keys()]):
        bulk_group = f['profiles'][str(lognum)]['bulk']
        model_number = f['profiles'][str(lognum)]['headers']['model_number'][()]
        print('reading profile {} / model number {}'.format(lognum, model_number))

        mesa_profiles = dict()
        for key in ['mass', 'kt', 'kmu', 'gradT_sub_grada', 'gradr_sub_grada', 'gradL_composition_term',\
                    'eps_nuc', 'mixing_type', 'radius', 'tau', 'zone',  'log_D_mix',]:
            mesa_profiles[key] = bulk_group[key][()]
        locals().update(mesa_profiles)

        kt = 10**kt
        kmu = 10**kmu

        radius *= 6.957e10 #convert radius to solar units

        R0 = compute_R0(gradT_sub_grada, gradr_sub_grada, gradL_composition_term, use_gradT = True)
        tau = compute_tau(kmu, kt)

        curr_mass = f['profiles'][str(lognum)]['headers']['star_mass'][()]
        mass_coord = mass/curr_mass

        max_burn = np.argmax(eps_nuc)
        good_mass = (mass_coord < mass_coord[max_burn]*1.1)
        if max_burn != len(eps_nuc) - 1:
            good_mass *= (mass_coord > mass_coord[max_burn+1])
        else:
            good_mass *= (mass_coord > mass_coord[max_burn])
        thermohaline = mixing_type==4
        good_bool = good_mass*thermohaline
        if np.sum(good_bool) >= 10:
            #remove outliers
            delta_m = np.gradient(mass[good_bool])
            good_bool[good_bool] *= np.abs(delta_m) < 3*np.abs(np.median(delta_m))

            #volume average of R0 and tau
            radius_th = radius[good_bool]
            mass_th = mass_coord[good_bool]
            xmin_th = mass_th.min()
            xmax_th = mass_th.max()
            ymin = np.min(np.log10(eps_nuc[good_mass])) - 3
            ymax = np.max(np.log10(eps_nuc[good_mass])) + 5
            radial_extent = radius_th.max() - radius_th.min()
            mass_extent = xmax_th - xmin_th
            drs.append(radial_extent)
            dms.append(mass_extent)


            #Only look at 1/3 of the zone, near the bottom
            good_bool *= (mass_coord > mass_th.min() + 0.1*mass_extent)*(mass_coord <= mass_th.min() + 0.433*mass_extent)
#            good_bool *= (radius > radius_th.min() + 0.25*radial_extent)*(radius <= radius_th.max() - 0.25*radial_extent)
            radius_th = radius[good_bool]
            mass_th = mass_coord[good_bool]

            R0_th = R0[good_bool]
            tau_th = tau[good_bool]
            Dmix = 10**(log_D_mix[good_bool])
            dV = 4*np.pi*radius_th**2
            thermohaline_r = (R0_th - 1)/(tau_th**(-1) - 1)
            avg_R0 = np.trapz(R0_th*dV, x=radius_th)/np.trapz(dV, x=radius_th)
            avg_tau = np.trapz(tau_th*dV, x=radius_th)/np.trapz(dV, x=radius_th)
            avg_th_r = np.trapz(thermohaline_r*dV, x=radius_th)/np.trapz(dV, x=radius_th)
            avg_Dmix = np.trapz(Dmix*dV, x=radius_th)/np.trapz(dV, x=radius_th)

            good_ages.append(f['profiles'][str(lognum)]['headers']['star_age'][()])
            good_log_g.append(log_gs[star_age == good_ages[-1]])
            good_models.append(model_number)
            good_R0s.append(avg_R0)
            good_taus.append(avg_tau)
            good_rs.append(avg_th_r)

            tau_mixs.append(drs[-1]**2/avg_Dmix)
            if len(good_point) < 21:
                good_point.append(False)
            else:
                diff = np.diff(dms[-21:])
                #Once zone size stops changing appreciably, mark measures as good.
                avg_diff = np.abs(np.mean(diff)/np.abs(dms).max()) 
                print('avg_diff: {:.3e}'.format(avg_diff))
                if avg_diff < 5e-3:
                    good_point.append(True)
                else:
                    good_point.append(False)
            log_str = 'age: {:.2e}, tau_mix: {:.2e}, \ndr = {:.3e} cm, dm = {:.2e}, good = {}'.format(good_ages[-1], tau_mixs[-1], radial_extent, mass_extent, good_point[-1])

            print(log_str)
            if PLOT_MOVIE:
                ax.set_prop_cycle('color', palettable.colorbrewer.qualitative.Dark2_8.mpl_colors)
                ax.fill_between((xmin_th, xmax_th), ymin-2, ymax+2, color='grey', alpha=0.25)
                #update for trimmed edges
                xmin_th = mass_th.min()
                xmax_th = mass_th.max()
                ax.fill_between((xmin_th, xmax_th), ymin-2, ymax+2, color='grey', alpha=0.25)
                ax.plot(mass_coord[good_mass], np.log10(R0)[good_mass], label='log10(R0)')
                ax.plot(mass_coord[good_mass], np.log10(1/tau)[good_mass], label='log10(1/tau)')
                ax.plot(mass_coord[good_mass], np.log10(eps_nuc)[good_mass], label='log10(eps_nuc)')
                ax.plot(mass_coord[good_mass], log_D_mix[good_mass], label='log_D_mix')

                plt.suptitle(log_str)

                axins.plot(mass_coord[good_bool], np.log10(R0_th), c=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors[0])
                axins.plot(mass_coord[good_bool], np.log10(1/tau_th), c=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors[1])
                axins.axhline(np.log10(avg_R0), c=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors[0])
                axins.axhline(np.log10(1/avg_tau), c=palettable.colorbrewer.qualitative.Dark2_8.mpl_colors[1])
                ax.legend(loc='lower right')
                ax.set_ylim(ymin, ymax)
#                ax.set_xscale('log')
                fig.savefig('figs_R0/R0_fig_{:06d}.png'.format(lognum))
                ax.clear()
                axins.clear()
        if np.sum(good_point) > 1000:
            break


good_models = np.array(good_models)
good_ages = np.array(good_ages)
good_R0s = np.array(good_R0s)
good_taus = np.array(good_taus) 
good_rs = np.array(good_rs)
good_log_g = np.array(good_log_g)


drs = np.array(drs)
dms = np.array(dms)
tau_mixs = np.array(tau_mixs)
good_point = np.array(good_point, dtype=bool)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
print(good_models, good_R0s)

ax.set_prop_cycle('color', palettable.colorbrewer.qualitative.Dark2_8.mpl_colors)
ax.plot(good_models, np.log10(np.array(good_R0s)), label='log10(R0)')
ax.plot(good_models, np.log10(1/np.array(good_taus)), label='log10(1/tau)')
ax.set_xlabel('model number')
ax.set_ylabel('various')
ax.legend()
fig.savefig('figs_R0/full_R0_vs_model.png', dpi=300)
ax.clear()

ax.plot(models, log_gs)
ax.set_xlabel('model number')
ax.set_ylabel('log_g')
fig.savefig('figs_R0/full_log_g_vs_model.png', dpi=300)
ax.clear()

r_post = (good_R0s - 1)/(good_taus**(-1) - 1)
ax.plot(good_models, good_rs, label='during', lw=2)
ax.plot(good_models, r_post, label='post')
ax.set_ylim(1e-6, 1)
ax.set_yscale('log')
ax.set_xlabel('model number')
ax.set_ylabel(r'$r = (R_0 - 1)/(\tau^{-1} - 1)$')
ax.legend()
fig.savefig('figs_R0/full_r_vs_model.png', dpi=300)
ax.clear()

ax.plot(good_log_g, good_rs, label='during', lw=2)
ax.plot(good_log_g, r_post, label='post')
ax.set_ylim(1e-6, 1)
ax.set_xlabel('log g')
ax.set_ylabel(r'$r = (R_0 - 1)/(\tau^{-1} - 1)$')
ax.set_yscale('log')
ax.legend()
fig.savefig('figs_R0/full_r_vs_log_g.png', dpi=300)
ax.clear()


with h5py.File('r_vs_time.h5', 'w') as f:
    f['models'] = good_models
    f['ages'] = good_ages
    f['r'] = good_rs
    f['R0'] = good_R0s
    f['tau'] = good_taus
    f['log_g'] = good_log_g

    f['dr'] = drs
    f['dm'] = dms
    f['tau_mix'] = tau_mixs
    f['good_point'] = good_point

