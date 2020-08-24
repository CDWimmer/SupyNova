from arnett_small_R import bolo_l
from nagy.nagy import kappa_nagy
from lmfit import Model, Parameters
from collections import Iterable
from numpy import loadtxt, log10, array, sqrt
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

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

data = "data/2009jf_bolom.txt"
data_skip = 2
phase, log_lum, uncert = loadtxt(data, unpack=True)

time = array([(t + abs(phase[0]))*day2sec for t in phase])[data_skip:]
log_lum = log_lum[data_skip:]
# time[0]=time[1]/2
# non-zero start or bad things tend happen
print("t[0] =", time[0])

lc_model = Model(fit_me)
print('parameter names: {}'.format(lc_model.param_names))
print('independent variables: {}'.format(lc_model.independent_vars))

# params = lc_model.make_params(m=6 * m_sun, v_sc=1e9, m_ni=0.7 * m_sun, shift=500)
params = Parameters()
params.add(name="m", value=1*m_sun, min=1*m_sun, max=20*m_sun)
params.add(name="v_sc", value=5e9, min=5e8, max=20e9)
params.add(name="m_ni", value=0.1*m_sun, min=0.001*m_sun, max=1*m_sun)
params.add(name="shift", value=500, min=1, max=3*day2sec)


# test = lc_model.eval(params=params, times=time)
# plt.plot(time/day2sec, log_lum, label="Data")
# plt.plot(time/day2sec, test, label="Initial Parameters")
# plt.xlabel("time (days)")
# plt.ylabel(r"log$_{10}$ luminosity")
# plt.legend()
# plt.show()

# Actually do a fit:
print("Fitting....")
try:
    result = lc_model.fit(log_lum, params, times=time)
except Exception as e:
    print(f"Something went wrong during fitting, perhaps parameters went mad: \033[91m{str(e)}")
    exit(1)
else:
    print("Done!")
    print(result.fit_report())  # gives things like chi square
    fig, ax = plt.subplots(1, 1)
    ax.semilogx(time/day2sec, log_lum, 'bo', label=f"Data for {data.split('/')[1].split('.')[0].split('_')[0]}")
    ax.plot(time/day2sec, result.init_fit, 'k--', label="Initial fit")
    ax.plot(time/day2sec, result.best_fit, 'r-', label="Best fit")
    ax.set_xlabel("time (days)")
    ax.set_ylabel(r"log$_{10}$ luminosity")
    ax.legend()
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    ax.xaxis.set_major_formatter(formatter)
    plt.show()

print(f"M = {result.params['m'].value/m_sun:.5f} solar masses\n"
      f"v_sc = {result.params['v_sc'].value/1e9:.5f} e9 cm/s\n"
      f"M_Ni = {result.params['m_ni'].value/m_sun:.5f} solar masses")

print(f"sqrt(2.53M_Ni/M) = {sqrt(2.53*result.params['m_ni'].value/result.params['m'].value)}")

# todo: Do a finer plot with the fitted parameters (e.g. more and evenly-spaced timesteps:
