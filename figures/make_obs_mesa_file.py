import glob
import numpy as np

base_obs_file = 'mixing_vs_r/scratchtableforEvan_figD.txt'
output_file = 'mixing_vs_r/observed_mixing_and_r.csv'

data = np.genfromtxt("./mixing_vs_r/scratchtableforEvan_figD.txt", skip_header=4)
csvfiles = glob.glob('mesa_spread/*.csv')

M = data[:,0]
feh = data[:,1]
mixing_apk = data[:,2]
mixing_shet = data[:,3]

rs = np.zeros((4, len(mixing_apk)))
R0s = np.zeros((4, len(mixing_apk)))
taus = np.zeros((4, len(mixing_apk)))


for fname in csvfiles:
    if 'Brown' in fname:
        i = 0
    elif 'coeff1.0e-01' in fname:
        i = 1
    elif 'coeff2.0e+00' in fname:
        i = 2
    elif 'coeff7.0e+02' in fname:
        i = 3


    data = np.genfromtxt(fname, delimiter=',', skip_header=1)
    mesa_M = data[:,0]
    mesa_FeH = data[:,1]
    mesa_r = data[:,3]
    mesa_R0 = data[:,4]
    mesa_tau = data[:,5]

    for j in range(len(M)):
        goodval = (mesa_M == M[j])*(mesa_FeH == feh[j])
        if np.sum(goodval) == 1:
            rs[i,j] = mesa_r[goodval]
            R0s[i,j] = mesa_R0[goodval]
            taus[i,j] = mesa_tau[goodval]
        else:
            rs[i,j] = np.nan
            R0s[i,j] = np.nan
            taus[i,j] = np.nan

data = []

for i in range(len(M)):
    this_data = [M[i], feh[i], mixing_apk[i], mixing_shet[i],]
    this_data += [rs[j,i] for j in range(4)]
    this_data += [R0s[j,i] for j in range(4)]
    this_data += [taus[j,i] for j in range(4)]
    data.append(this_data)
    print(data[-1])
header = "{:>21s}".format("M") + (15*"{:>23s}").format("Fe/H", "DeltaCN_corr_APOKASC", "DeltaCN_corr_Shetrone", \
        "r_brown_1", "r_kip_0.1", "r_kip_2", "r_kip_700",\
        "R0_brown_1", "R0_kip_0.1", "R0_kip_2", "R0_kip_700",\
        "tau_brown_1", "tau_kip_0.1", "tau_kip_2", "tau_kip_700")

np.savetxt(output_file, data, '%22.10f', delimiter=',', header=header)

print(rs)
