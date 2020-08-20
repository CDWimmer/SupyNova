from arnett_small_R import bolo_l
from nagy.nagy import kappa_nagy
from lmfit import Model, Parameters
from collections import Iterable
from numpy import loadtxt, log10, array
from matplotlib import pyplot as plt

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams


def fit_me(times: Iterable, m, v_sc, m_ni, shift):
    """
    Return bolometric luminosity values for a given set of times.
    Uses the kappa_nagy function for kappa calculation.
    """
    time_shifted = [t+shift for t in times]
    results = []
    for t in time_shifted:
        results.append(log10(bolo_l(t, m, v_sc, m_ni, kappa=kappa_nagy)))
    return results


data_skip = 1
phase, log_lum, uncert = loadtxt("data/2006aj_nature.txt", unpack=True)

time = array([(t + abs(phase[0]))*day2sec for t in phase])[data_skip:]
log_lum = log_lum[data_skip:]
# time[0]=time[1]/2
# non-zero start or bad things tend happen


lc_model = Model(fit_me)
print('parameter names: {}'.format(lc_model.param_names))
print('independent variables: {}'.format(lc_model.independent_vars))

# params = lc_model.make_params(m=6 * m_sun, v_sc=1e9, m_ni=0.7 * m_sun, shift=500)
params = Parameters()
params.add(name="m", value=6*m_sun, min=0.75*m_sun, max=20*m_sun)
params.add(name="v_sc", value=5e9, min=1e8, max=20e9)
params.add(name="m_ni", value=0.7*m_sun, min=0.1*m_sun, max=2*m_sun)
params.add(name="shift", value=500, min=0, max=3*day2sec)


test = lc_model.eval(params=params, times=time)
plt.plot(time/day2sec, log_lum, label="Data")
plt.plot(time/day2sec, test, label="Initial Parameters")
plt.xlabel("time (days)")
plt.ylabel(r"log$_{10}$ luminosity")
plt.legend()
plt.show()

# Actually do a fit:
print("Fitting....")
try:
    result = lc_model.fit(log_lum, params, times=time)
except Exception as e:
    print(f"Something went wrong during fitting, perhaps parameters went mad: \033[91m{str(e)}")
else:
    print("Done!")
    print(result.fit_report())  # gives things like chi square
    plt.semilogx(time/day2sec, log_lum, 'bo')
    plt.plot(time/day2sec, result.init_fit, 'k--', label="initial fit")
    plt.plot(time/day2sec, result.best_fit, 'r-', label="best fit")
    plt.xlabel("time (days)")
    plt.ylabel(r"log$_{10}$ luminosity")
    plt.legend()
    plt.show()


# todo: Do a finer plot with the fitted parameters (e.g. more and evenly-spaced timesteps:
