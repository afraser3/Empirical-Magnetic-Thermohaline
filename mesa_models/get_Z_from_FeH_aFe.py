#!/usr/bin/env python3
# Written by Tim White, Meridith Joyce, Evan Anders
# Uses code from Aaron Dotter's repo: https://github.com/aarondotter/initial_xa_calculator/blob/master/src/initial_xa_calculator.f90#L303
import numpy as np
import pandas as pd
import argparse

dYdZ = 1.3426
big_bang_Y = 0.2485
z_to_y = lambda zval: big_bang_Y + dYdZ * zval

def pow10(A):
    return 10**(A)

#########################
def numfrac_justZ(feh,alpha,abundances):
    table = pd.read_fwf('GS98abundances.txt')

    setup(feh,alpha,table)

    b = pow10(table['A'])
    errb = b*np.log(10)*table['errA']
    sum0 = np.sum(b)
    errsum0 = np.sqrt(np.sum((errb)**2))

    b1 = np.copy(b)
    b1[0] = 0.0
    b1[1] = 0.0

    sum1 = np.sum(b1)

    b = b/sum0
    sumam = np.sum(b*table['awt'])
    b1 = b1/sum1
    errb1 = errb/sum1
    errb = errb/sum0

    a = b*table['awt']/sumam
    erra = errb*table['awt']/sumam
    
    sumz = np.sum(a[2:])
    errsumz = np.sqrt(np.sum(erra[2:]**2))  

    return sumz, errsumz, sumz/a[0], errsumz/a[0]

###############

def setup(feh,alpha,table):
    for index, row in table.iterrows():
        if row['element'] == 'C':
            table.at[index,'A'] += feh
        if row['element'] == 'N':
            table.at[index,'A'] += feh
        if row['element'] == 'O':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'Ne':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'Na':
            table.at[index,'A'] += feh
        if row['element'] == 'Mg':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'Al':
            table.at[index,'A'] += feh
        if row['element'] == 'Si':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'P':
            table.at[index,'A'] += feh
        if row['element'] == 'S':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'Cl':
            table.at[index,'A'] += feh
        if row['element'] == 'Ar':
            table.at[index,'A'] += feh
        if row['element'] == 'Ca':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'Ti':
            table.at[index,'A'] += feh
            table.at[index,'A'] += alpha
        if row['element'] == 'Cr':
            table.at[index,'A'] += feh
        if row['element'] == 'Mn':
            table.at[index,'A'] += feh
            table.at[index,'A'] -= alpha
        if row['element'] == 'Fe':
            table.at[index,'A'] += feh
        if row['element'] == 'Ni':
            table.at[index,'A'] += feh

if __name__ == '__main__':
    abundances = 'GS98' 

    feh_array = np.array([0.4, 0.2, 0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -1.2, -1.4])
    aFe_array = [0.0] 

    xvals = []
    yvals = []
    zvals = []

    print((10*"{:<10s}").format('FeH,', 'alpha/Fe,', 'sumZ,', 'err_sumZ,', 'Z/X,', 'Z/X_err,', 'X', 'Y', 'Z', 'Y_func'))
    for feh in feh_array:
        for alpha in aFe_array:
            sumZ, err_sumZ, ZX, err_ZX = numfrac_justZ(feh, alpha, abundances)
            new_X = (1 - big_bang_Y)/(1 + ZX*(1 + dYdZ))
            new_Z = ZX * new_X
            new_Y = 1 - (new_X + new_Z)
            zvals.append(new_Z)
            yvals.append(new_Y)
            xvals.append(new_X)
            Y_round = z_to_y(new_Z)
            print((10*'{:<10.4f}').format(feh, alpha, sumZ, err_sumZ, ZX, err_ZX, new_X, new_Y, new_Z, Y_round))

    z_over_x = np.log10(np.array(zvals)/np.array(xvals))
    z_over_x -= z_over_x[feh_array == 0]
    rows = [[feh, x, y, z, log_z_over_x] for feh, x, y, z, log_z_over_x in zip(feh_array, xvals, yvals, zvals, z_over_x)]
    header = "{:>9s}".format("Fe/H") + (4*"{:>12s}").format("X", "Y", "Z", "log10(Z/X)")
    np.savetxt('Z_feH_table.dat', rows, '%11.8f', header=header)
