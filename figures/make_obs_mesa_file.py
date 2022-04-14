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

    for j in range(len(M)):
        goodval = (mesa_M == M[j])*(mesa_FeH == feh[j])
        if np.sum(goodval) == 1:
            rs[i,j] = mesa_r[goodval]
        else:
            rs[i,j] = np.nan

data = []

for i in range(len(M)):
    data.append((M[i], feh[i], mixing_apk[i], mixing_shet[i], rs[0,i], rs[1,i], rs[2,i], rs[3,i]))
    print(data[-1])
header = "{:>18s}".format("M") + (7*"{:>21s}").format("Fe/H", "DeltaCN_corr_APOKASC", "DeltaCN_corr_Shetrone", "r_brown_1", "r_kip_0.1", "r_kip_2", "r_kip_700")

np.savetxt(output_file, data, '%20.8f', delimiter=',', header=header)

print(rs)
