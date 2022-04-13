import numpy as np
from scipy import optimize as opt
import matplotlib
import matplotlib.pyplot as plt
from palettable.colorbrewer.qualitative import Dark2_4
plt.style.use('../apj.mplstyle')

def rfromR(R0, tau):
    """
    Calculates the "reduced density ratio" r from the actual density ratio R0.
    See Brown et al. 2013 Eq. 11
    """
    return (R0 - 1.0) / (-1.0 + 1.0 / tau)


def Rfromr(r, tau):
    return 1.0 + r * (-1.0 + 1.0 / tau)


def lamguess(pr, tau, R0):
    """
    Guesses lamda_FGM (growth rate of the fastest-growing mode) according to Brown et al. 2013 Appendix B
    """
    r = rfromR(R0, tau)
    if r < tau:
        return np.sqrt(pr) - pr * np.sqrt(1.0 + tau / pr)
    else:
        if r > 0.5:
            return 2.0 * pr * (tau / pr) * ((1.0 / 3.0) * (1.0 - r)) ** (3.0 / 2.0) / (
                    1.0 - (1.0 - r) * (1.0 + tau / pr) / 3.0)
        else:
            return np.sqrt(pr * tau / r) - pr * np.sqrt(1 + tau / pr)


def k2guess(pr, tau, R0):
    """
    Guesses the wavenumber (squared) of the fastest-growing mode, see the docstring for lamguess above
    """
    r = rfromR(R0, tau)
    if r < tau:
        return (1.0 + tau / pr) ** (-0.5) - np.sqrt(pr) * (1.0 + (tau / pr) * (1.0 + tau / pr) ** (-2.0))
    else:
        if r > 0.5:
            return np.sqrt((1.0 - r) / 3.0)
        else:
            return np.sqrt((1.0 + tau / pr) ** (-0.5) - 2.0 * np.sqrt(r * tau / pr) * (1.0 + tau / pr) ** (-5.0 / 2.0))


def eq1(lam, k2, pr, tau, R0):
    """
    This function, as well as the next three functions (eq2, fun, and jac) are needed by gaml2max for calculating the
    fastest-growing mode's growth rate and wavenumber according to Sec. 4.1 of Brown et al. 2013
    """
    b2 = k2 * (1.0 + pr + tau)
    b1 = k2 ** 2.0 * (tau * pr + pr + tau) + pr * (1.0 - 1.0 / R0)
    b0 = k2 ** 3.0 * tau * pr + k2 * pr * (tau - 1.0 / R0)
    return lam ** 3.0 + b2 * lam ** 2.0 + b1 * lam + b0


def eq2(lam, k2, pr, tau, R0):
    c2 = 1.0 + pr + tau
    c1 = 2.0 * k2 * (tau * pr + tau + pr)
    c0 = 3.0 * k2 ** 2.0 * tau * pr + pr * (tau - 1.0 / R0)
    return c2 * lam ** 2.0 + c1 * lam + c0


def fun(x, pr, tau, R0, passk1=False):  # returns f(x) where f = [eq1, eq2] and x = [lam, k2]
    if passk1:  # if x[1] is k instead of k^2
        return [eq1(x[0], x[1] ** 2.0, pr, tau, R0), eq2(x[0], x[1] ** 2.0, pr, tau, R0)]
    else:
        return [eq1(x[0], x[1], pr, tau, R0), eq2(x[0], x[1], pr, tau, R0)]


def jac(x, pr, tau, R0, passk1=False):  # jacobian of fun(x)
    lam = x[0]
    if passk1:  # is x[1] k or k^2?
        k2 = x[1] ** 2.0
    else:
        k2 = x[1]
    b2 = k2 * (1.0 + pr + tau)
    db2dk2 = 1.0 + pr + tau  # derivative of b2 wrt k2
    b1 = k2 ** 2.0 * (tau * pr + pr + tau) + pr * (1.0 - 1.0 / R0)
    db1dk2 = 2.0 * k2 * (tau * pr + pr + tau)
    b0 = k2 ** 3.0 * tau * pr + k2 * pr * (tau - 1.0 / R0)
    db0dk2 = 3.0 * k2 ** 2.0 * tau * pr + pr * (tau - 1.0 / R0)

    j11 = 3.0 * lam ** 2.0 + 2.0 * b2 * lam + b1  # d(eq1)/dlam
    j12 = lam ** 2.0 * db2dk2 + lam * db1dk2 + db0dk2  # d(eq1)/dk2
    if passk1:
        j12 = j12 * 2.0 * x[1]  # d(eq1)/dk = d(eq1)/dk2 * dk2/dk

    c2 = 1.0 + pr + tau
    c1 = 2.0 * k2 * (tau * pr + tau + pr)
    dc1dk2 = c1 / k2
    c0 = 3.0 * k2 ** 2.0 * tau * pr + pr * (tau - 1.0 / R0)
    dc0dk2 = 6.0 * k2 * tau * pr

    j21 = 2.0 * c2 * lam + c1
    j22 = lam * dc1dk2 + dc0dk2
    if passk1:
        j22 = j12 * 2.0 * x[1]
    return [[j11, j12], [j21, j22]]


def gaml2max(pr, tau, R0):
    """
    Uses scipy.optimize.root with the above functions to find lambda_FGM and l^2_FGM.
    See section 4.1 of Brown et al. 2013
    """
    sol = opt.root(fun, [lamguess(pr, tau, R0), k2guess(pr, tau, R0)], args=(pr, tau, R0), jac=jac, method='hybr')
    x = sol.x
    if sol.x[1] < 0:  # if a negative k^2 is returned, then try again but solve for k instead of k^2
        sol = opt.root(fun, [lamguess(pr, tau, R0), np.sqrt(k2guess(pr, tau, R0))], args=(pr, tau, R0, True), jac=jac,
                       method='hybr')
        test = fun(sol.x, pr, tau, R0, True)
        if np.allclose(test, np.zeros_like(test)) == False:
            raise ValueError("gaml2max is broken!")
        x = sol.x
        x[1] = x[1] ** 2.0  # whatever calls gaml2max expects k^2, not k
    # sol = opt.root(fun, [lamguess(pr, tau, R0), 4.0*k2guess(pr, tau, R0)], args=(pr, tau, R0), jac=jac, method='hybr')
    # if sol.x[1]<0:
    # raise ValueError("gaml2max settled on a negative l2!")
    # return sol.x
    return x


def HG19_eq32(w, pr, tau, R0, HB, CH=1.66):
    """
    Simply evaluates Eq. (32) in Harrington & Garaud 2019.
    Specifically, evaluates LHS - RHS (so it should evaluate to zero if w is the solution)

    Parameters
    ----------
    w: w_f in the equation
    pr: Prandtl number, used to calculate lambda_f and l_f
    tau: ratio of diffusivities, used to calculate lambda_f and l_f
    R0: density ratio, used to calculate lambda_f and l_f
    HB: Lorentz force coefficient
    CH: fitting parameter

    Returns
    -------
    LHS - RHS of Eq. (32)
    """

    lamhat, l2hat = gaml2max(pr, tau, R0)
    lhat = np.sqrt(l2hat)
    LHS = 0.5 * w**2.0 - HB
    RHS1 = (CH*lamhat/(0.42*lhat))**1.5
    RHS2 = np.sqrt(w)
    #print(w)
    return LHS - (RHS1*RHS2)
    # return 0.5 * w**2.0 - HB - (CH*lamhat/(0.42*lhat))**1.5 * np.sqrt(w)


def dEQ32dw(w, pr, tau, R0, HB, CH=1.66):
    """
    Derivative with respect to w of the function HG19_eq32 above. For list of inputs, see HG19_eq32.
    Where HG19_eq32 evaluates F(w) = LHS - RHS (LHS and RHS of Eq. 32), so that
    F(w) = 0 for the solution w, this function returns dF/dw for use in root-finding algorithm
    """

    lamhat, l2hat = gaml2max(pr, tau, R0)
    lhat = np.sqrt(l2hat)
    return w - 0.5 * (CH*lamhat/(0.42*lhat))**1.5 / np.sqrt(w)


def w_f_HG19(pr, tau, R0, HB, CH=1.66):
    """
    Uses a root-finding algorithm to solve EQ32 in HG19 for w

    Parameters
    ----------
    pr: Prandtl number, used to calculate lambda_f and l_f
    tau: ratio of diffusivities, used to calculate lambda_f and l_f
    R0: density ratio, used to calculate lambda_f and l_f
    HB: Lorentz force coefficient
    CH: fitting parameter

    Returns
    -------
    w: the shear velocity that solves Eq. 32
    """

    lamhat, l2hat = gaml2max(pr, tau, R0)
    lhat = np.sqrt(l2hat)
    # the following is a terrible initial guess
    # w0 = np.sqrt(2.0*HB)
    w0 = max(np.sqrt(2.0*HB), 2.0 * np.pi * lamhat/lhat)
    # if HB > 1.0:
        # w0 = np.sqrt(2.0*HB)
    # else:
        # w0 = 2.0**(2.0/3.0) * CH*lamhat/(0.42*lhat)
    result = opt.root_scalar(HG19_eq32, args=(pr, tau, R0, HB, CH), x0=w0, fprime=dEQ32dw)
    root = result.root
    if root > 0:
        return result
    else:
        w1 = 10.0 * CH**1.5 * (np.sqrt(2.0*HB) + 2.0 * np.pi * lamhat/lhat)
        try:
            result = opt.root_scalar(HG19_eq32, args=(pr, tau, R0, HB, CH),
                                     bracket=[0.0, w1])
        except ValueError:
            print("w1 = ", w1)
            print("CH = ", CH)
            raise
        return result


def NuC_from_w(tau, w, lamhat, l2hat, KB=1.24):
    return 1 + KB * w ** 2.0 / (tau * (lamhat + tau * l2hat))


def NuC_Brown_model(pr, tau, R0):
    C_brown = 7.0096  # the C in Brown et al. 2013 Eq. 31 (they did a chi-squared it and got this value)
    lamhat, l2hat = gaml2max(pr, tau, R0)
    return 1 + C_brown**2 * lamhat**2 / (tau * l2hat * (lamhat + tau*l2hat))


def NuC_Traxler_model(pr, tau, R0):
    """
    Evaluates Brown et al. 2013 Eq. 17
    """
    r = rfromR(R0, tau)
    a_Traxler = 101
    b_Traxler = 3.6
    c_Traxler = 1.1
    g_Traxler = a_Traxler * np.exp(-b_Traxler*r) * (1 - r)**c_Traxler
    return 1 + np.sqrt(pr/tau)*g_Traxler


def NuC_Kippenhahn_model(tau, R0, alpha_th):
    return 1.0 + (3/2)*alpha_th/(tau*R0)




data = np.genfromtxt("../mixing_vs_r/scratchtableforEvan_figD.txt", skip_header=4)

M = data[:,0]
feh = data[:,1]
mixing_apk = data[:,2]
mixing_shet = data[:,3]
r = r_brown = data[:,4]

norm = matplotlib.colors.Normalize(vmin=feh.min(), vmax=feh.max())
sm = plt.cm.ScalarMappable(cmap='plasma', norm=norm)
sm.set_array([])

fig = plt.figure(figsize=(6.5, 2))
ax1 = fig.add_axes((0.00, 0.00, 0.45, 1.00))
ax2 = fig.add_axes((0.55, 0.00, 0.45, 0.92))
cax = fig.add_axes((0.55, 0.95, 0.45, 0.04)) 


#fig 1 reproduction
Pr = 1e-6
tau = 1e-6
rs = np.logspace(-5, 0)
R0s = Rfromr(rs, tau)

NuC_Brown = np.array([NuC_Brown_model(Pr, tau, r0) for r0 in R0s])
NuC_Kippenhahn_Cp2 = NuC_Kippenhahn_model(tau, R0s, 2.0)

HB2 = 1e-6
wfs_HG19_HB2 = [w_f_HG19(Pr, tau, r0, HB2).root for r0 in R0s]
NuCs_HG19_HB2 = np.zeros_like(wfs_HG19_HB2)
for i in range(len(wfs_HG19_HB2)):
    lamhat, l2hat = gaml2max(Pr, tau, R0s[i])
    NuCs_HG19_HB2[i] = NuC_from_w(tau, wfs_HG19_HB2[i], lamhat, l2hat)

#Plot nice background shading
r_models_min = 1e-4
r_models_max = 1e-3
N = 1000
color_spots = np.logspace(np.log10(r_models_min), np.log10(r_models_max), N)
y_spots = np.logspace(-10, 10, 2)

yy, xx = np.meshgrid(y_spots, color_spots)
color = np.log10(xx)
ax1.pcolormesh(xx, yy, color, vmin=np.log10(r_models_min), vmax=np.log10(r_models_max), cmap='viridis', alpha=0.5)

#Plot lines
ax1.loglog(rs, NuC_Brown-1,  c=Dark2_4.mpl_colors[1], label='Brown 19', lw=2)
ax1.loglog(rs, NuC_Kippenhahn_Cp2-1, c=Dark2_4.mpl_colors[0], label=r'Kippenhahn ($\alpha_{\rm{th}} = 2$)', lw=2)
ax1.loglog(rs, NuCs_HG19_HB2-1, c=Dark2_4.mpl_colors[2], label=r'HG19 ($H_B = 10^{-6}$)', lw=2, zorder=1)
ax1.set_xlim(1e-5, 1)
ax1.set_xlabel(r'$r$')
ax1.set_ylabel(r'$D_{\rm{th}}/\kappa_\mu$')
ax1.set_ylim(1e-1, 1e7)
ax1.legend(fontsize=8, frameon=True)


#fig 4 reproduction
kwargs = {'edgecolor' : 'k'}
for m in range(len(M)):
    ax2.scatter(r[m], mixing_shet[m], marker='s', color=sm.to_rgba(feh[m]), label='Shetrone', **kwargs)
    ax2.scatter(r[m], mixing_apk[m], marker='o', color=sm.to_rgba(feh[m]), label='APOGEE', **kwargs)
    if m == 0:
        ax2.legend(loc='center right', fontsize=8)
ax2.text(0.98, 0.98, 'Brown', ha='right', va='top', transform=ax2.transAxes, fontsize=9)
ax2.set_xlim(1e-5, 1)
ax2.set_ylim(-0.05, 0.7)
ax2.set_xscale('log')

cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=matplotlib.cm.get_cmap('plasma'), norm=norm, orientation='horizontal', ticklocation='top')
cbar.set_label(r'[Fe/H]')

ax1.set_xlabel(r'$r$')
ax2.set_xlabel(r'$r$')
ax2.set_ylabel(r'$\Delta$ [C/N]')


fig.savefig('punchline.png', dpi=300, bbox_inches='tight')
fig.savefig('punchline.pdf', dpi=300, bbox_inches='tight')



