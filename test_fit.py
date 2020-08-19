from arnett_small_R import bolo_l
from nagy.nagy import kappa_nagy
from lmfit import Model
from collections import Iterable
from numpy import loadtxt, log10, array
from matplotlib import pyplot as plt

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams


def fit_me(times: Iterable, m, v_sc, m_ni):
    """
    Return bolometric luminosity values for a given set of times.
    Uses the kappa_nagy function for kappa calculation.
    """
    results = []
    for t in times:
        results.append(log10(bolo_l(t, m, v_sc, m_ni, kappa=kappa_nagy)))
    return results


phase, log_lum, uncert = loadtxt("data/2011fe_bolom.txt", unpack=True)

time = array([(t + abs(phase[0]))*day2sec for t in phase])[1:]
log_lum = log_lum[1:]
#time[0]=time[1]/2
# non-zero start or bad things tend happen


lc_model = Model(fit_me)
print('parameter names: {}'.format(lc_model.param_names))
print('independent variables: {}'.format(lc_model.independent_vars))
params = lc_model.make_params(m=6 * m_sun, v_sc=1e9, m_ni=0.7 * m_sun)  # about 10% nickel, probably a lot but we'll see

test = lc_model.eval(params=params, times=time)
plt.plot(time/day2sec, log_lum, label="Data")
plt.plot(time/day2sec, test, label="Initial Parameters")
plt.xlabel("time (days)")
plt.ylabel(r"log$_{10}$ luminosity")
plt.legend()
plt.show()

# Actually do a fit:
print("Fitting....")
result = lc_model.fit(log_lum, params, times=time)
print("Done!")
print(result.fit_report())  # gives things like chi square
plt.semilogx(time/day2sec, log_lum, 'bo')
plt.plot(time/day2sec, result.init_fit, 'k--', label="initial fit")
plt.plot(time/day2sec, result.best_fit, 'r-', label="best fit")
plt.xlabel("time (days)")
plt.ylabel(r"log$_{10}$ luminosity")
plt.legend()
plt.show()


