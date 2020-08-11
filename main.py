import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from lmfit import Model


# common constants
c = 2.9979e+10 # cm/s
day = 86400  # s
M_sun = 1.989e+33  # solar masses to grams.
kappa_gamma = 0.03  # opacity of gamma rays
beta = 13.7  # constant used by Arnett

tau_56ni = 6.077 * day / np.log(2)  # decay time nickel 56
tau_56co = 77.27 * day / np.log(2)  # decay time cobalt 56



# Set up timescale parameter function
NAGY_KAPPA15 = "../data/15M_sun.csv"
NAGY_KAPPA30 = "../data/30M_sun.csv"
t_k, kappa_t = np.loadtxt(NAGY_KAPPA15, unpack=True, delimiter=", ")

low_t = min(t_k)
high_t = max(t_k)
k_inter = interp1d(t_k, kappa_t)  # linear interpolation - could try cubic or something?


def ts_nagy(t, vars):
    """Diffusion Timescale Parameter function"""
    E_k = 3 / 10 * vars["M_ej"] * (vars["V_ej"] ** 2)
    if t <= low_t:
        _kappa = kappa_t[0]  # first val
    elif t >= high_t:
        _kappa = kappa_t[-1]  # final val
    else:
        _kappa = k_inter(t)
    return 1.05 / np.sqrt(beta * c) * np.sqrt(_kappa) * vars["M_ej"] ** (3 / 4) * E_k ** (-1 / 4)


# Deposition function/s:

def deposition(t, vars):
    """Deposition function for nickel, Arnett eqn. 50/51
        Requires vars: timescale_parameters, V_ej"""
    R = vars["V_ej"] * ts_nagy(t, vars)
    # R = V_ej*timescale_parameter(t, vars) # should this take the timescale parameter or just time?
    rho = vars["M_ej"] / (4 / 3 * np.pi * R ** 3)
    tau_56co_gamma = kappa_gamma * rho * R  # diffusion timescale of photons of decay from ni->co.
    G = tau_56co_gamma / (tau_56co_gamma + 1.6)
    D = G*(1+2*G*(1-G)*(1-0.75*G))
    return D


# def d_cobalt56(t, vars):
#     """Requires vars: V_ej, M_ej"""
#     if not {"V_ej", "M_ej"}.issubset(vars):
#         raise ValueError("Missing required args.")
#     """Deposition function for cobalt"""
#     res = 0.966 * d_nickel56(t, vars) + 0.034
#     return res


# Luminosity sub-functions:

def l_nickel56(t, vars):
    # Arnett Eq 31:
    # MATLAB XXXXXXXXXX
    res = deposition(t, vars) * 2 * t * np.e ** ((-t * ts_nagy(t, vars) / tau_56ni) + (t ** 2))
    return res


def l_cobalt56(t, vars):

    l_cobalt = (0.966 * deposition(t, vars) + 0.034) * 2 * t * (np.e ** (-t * ts_nagy(t, vars) / tau_56co) - np.e ** (-t * ts_nagy(t, vars) / tau_56ni)) / (1 - (tau_56ni / tau_56co) * np.e ** (t ** 2))
    return l_cobalt


def luminosity(t, M_ej, V_ej, preval, ni56_multiplier):
    vars = {"M_ej": M_ej,
            "V_ej": V_ej,
            #"preval": preval,
            "M_56Ni": M_ej * ni56_multiplier}
    val1 = np.e ** (-(t / ts_nagy(t, vars)) ** 2)  # common to both functions
    integral1 = integrate.quad(l_nickel56,
                               a=0.000001,
                               b=t / ts_nagy(t, vars),
                               args=vars)[0]
    integral2 = integrate.quad(l_cobalt56,
                               a=0.000001,
                               b=t / ts_nagy(t, vars),
                               args=vars)[0]

    res = np.log10(
        val1 * preval * vars["M_56Ni"] * (integral1 + integral2)
    )
    return res

y = []
for i in range(3000):
    print(i)
    y.append(luminosity(i, M_ej=15 * M_sun, V_ej=1e9, preval=1e25, ni56_multiplier=0.01))

plt.plot(y)
plt.show()

