import numpy as np
from scipy.integrate import quad
from nagy import nagy

# constants
tau_ni = 7.605e5  # s | nickel-56 decay time
tau_co = 9.822e6  # s | cobalt-56 decay time

e_ni = 4.78e10  # ergs /g /s  | energy release of nickel decay
e_co = 2.56e8  # ergs /g /s   | energy release of cobalt decay (purely from positrons)

beta = 13.7  # approximately this for a variety of cases(A80) so will use throughout this project
kappa_gamma = 0.03  # opacity of gamma rays

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams
c = 2.99792458e10  # cm /s    | speed of light in CGS

# model parameters:
# t - time
# m - total ejected mass
# v_sc - velocity scale factor
# For large initial radius, R?


# functions
# nickel lambda
# cobalt lambda
# nickel bolometric luminosity - do we use the same for both?
# kappa - the opacity, import from nagy.py?

def kappa_const(t):
    """placeholder"""
    return 0.08


def l_ni(t, m, v_sc, kappa: callable = kappa_const):
    """Arnett eqn. 31 | Lambda(x,y)"""
    tau_m = np.sqrt((2 * kappa(t) * m) / (beta * c * v_sc))
    x = t / tau_m
    y = tau_m / (2 * tau_ni)

    def l_ni_integral(z):
        """z is the variable of integration"""
        res = np.exp((-2 * z * y) + pow(z, 2)) * 2 * z
        return res

    result = np.exp(-pow(x, 2)) * quad(l_ni_integral, 0, x)[0]
    print(f"x={x:.5f}\ty={y:.3f}\txy={(x * y):.3f}\tlambda_ni(x,y)={result:.3f}")
    return x * y, result


def l_co(t, m, v_sc, kappa: callable = kappa_const):
    """Arnet eqn. 38"""
    tau_m = np.sqrt((2 * kappa(t) * m) / (beta * c * v_sc))
    x = t / tau_m
    y = tau_m / (2 * tau_ni)
    yp = tau_m / (2 * tau_co)

    def l_co_integral(z):
        # MATLAB: 2*z.*(exp(-2*z*yp) - exp(-2*z*y))/(1 - tau_ni/tau_co).*exp(z.^2)).
        res = 2 * z * (np.exp(-2 * z * yp) - np.exp(-2 * z * y)) / (1 - tau_ni/tau_co) * np.exp(pow(z, 2))
        return res

    result = np.exp(-pow(x, 2)) * quad(l_co_integral, 0, x)[0]
    return result


def bolo_l(t, m, v_sc, m_ni, kappa: callable = kappa_const):
    #print(f"inputs: t={t}, m={m}, v_sc={v_sc}, m_ni={m_ni}")

    def deposition():
        tau_m = np.sqrt((2 * kappa(t) * m) / (beta * c * v_sc))
        R = v_sc * tau_m  # is this the correct V and tau?
        rho = m / (4 / 3 * np.pi * R ** 3)
        tau_56co_gamma = kappa_gamma * rho * R  # difussion timescale of photons of decay from ni->co. Is this named correctly?
        G = tau_56co_gamma / (tau_56co_gamma + 1.6)
        D = G * (1 + 2 * G * (1 - G) * (1 - 0.75 * G))
        return D
    # evaluate (D_ni * e_ni * m_ni * lambda_ni) + (D_co * e_co * m_ni * lambda_co)
    return (deposition() * e_ni * m_ni * l_ni(t, m, v_sc, kappa=kappa)[1]) + ((0.966 * deposition() + 0.034) * e_co *m_ni * l_co(t, m, v_sc, kappa=kappa))


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    """Do some test values"""
    skip_days = 1  # Days to skip at start. Shock breakout dominates initially
    end_day = 200  # day on which to end simulation
    # test values:
    m = 1.45 * m_sun
    v_sc = 1.2e9
    m_ni = 0.7 * m_sun

    times = np.linspace(skip_days * day2sec, end_day * day2sec, 200)

    # constant kappa
    lni_res = []
    lco_res = []
    xy_res = []
    lum_res = []

    for t in times:
        out = l_ni(t, m=m, v_sc=v_sc)
        xy_res.append(out[0])
        lni_res.append(out[1])
        lco_res.append(l_co(t, m=m, v_sc=v_sc))
    for t in times:
        lum_res.append(bolo_l(t, m=m, v_sc=v_sc, m_ni=0.7))

    # variable kappa
    lni_res_var = []
    lco_res_var = []
    xy_res_var = []
    lum_res_var = []

    for t in times:
        out = l_ni(t, m=m, v_sc=v_sc, kappa=nagy.kappa_nagy)
        xy_res_var.append(out[0])
        lni_res_var.append(out[1])
        lco_res_var.append(l_co(t, m=m, v_sc=v_sc, kappa=nagy.kappa_nagy))
    for t in times:
        lum_res_var.append(bolo_l(t, m=m, v_sc=v_sc, m_ni=0.7, kappa=nagy.kappa_nagy))

    # plotting
    print(f"Using test values of m=1.45*m_sun g, v_sc=1.2e9 cm/s ")
    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches(12, 6.4)

    ax[0].loglog(times / day2sec, lni_res, label=fr"$\Lambda_{{ni}}$ Constant kappa ($\kappa = {kappa_const(0)}$)")
    ax[0].loglog(times / day2sec, lni_res_var, label=r"$\Lambda_{{ni}}$ Variable kappa")
    ax[0].loglog(times / day2sec, lco_res, label=fr"$\Lambda_{{co}}$ Constant kappa ($\kappa = {kappa_const(0)}$)")
    ax[0].loglog(times / day2sec, lco_res_var, label=r"$\Lambda_{{co}}$ Variable kappa")
    ax[0].set_xlabel(r"time (days)")
    ax[0].set_ylabel(r"$\Lambda \left(x, y\right)$")

    ax[0].legend()

    ax[1].semilogy(times / day2sec, lum_res, label="Constant kappa")
    ax[1].semilogy(times / day2sec, lum_res_var, label="Variable kappa")
    ax[1].set_xlabel(r"time (days)")
    ax[1].set_ylabel(r"Luminosity (ergs s$^{-1}$)")

    ax[1].legend()

    plt.show()

# todo: Could move tau_m, x, y into bolo_l instead of recalculating in both cobalt and nickel parts.
