from arnett_small_R import bolo_l
from nagy.nagy import kappa_nagy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams

# test values:
m_default = 6 * m_sun
m_array = [1.5 * m_sun, 3 * m_sun, 4.5 * m_sun, 6 * m_sun, 7.5 * m_sun]
v_sc_default = 1e9
v_sc_array = [0.5e9, 1e9, 1.5e9, 2e9]
m_ni_default = 0.7 * m_sun
m_ni_array = [0.1 * m_sun, 0.25 * m_sun, 0.5 * m_sun, 0.7 * m_sun, 1 * m_sun]

skip_days = 1  # Days to skip at start. (shock breakout makes it inacurate)
end_day = 200  # day on which to end simulation
times = np.linspace(skip_days * day2sec, end_day * day2sec, 300)

m_results = dict()
v_sc_results = dict()
m_ni_results = dict()

for m in m_array:
    res = []
    for t in times:
        res.append(bolo_l(t, m, v_sc_default, m_ni_default, kappa=kappa_nagy))
    m_results[str(m)] = res

for v_sc in v_sc_array:
    res = []
    for t in times:
        res.append(bolo_l(t, m_default, v_sc, m_ni_default, kappa=kappa_nagy))
    v_sc_results[str(v_sc)] = res

for m_ni in m_ni_array:
    res = []
    for t in times:
        res.append(bolo_l(t, m_default, v_sc_default, m_ni, kappa=kappa_nagy))
    m_ni_results[str(m_ni)] = res


# prepare figures/axes

fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey="all")
fig.set_size_inches(10, 5)
ax1.set_xticklabels([])
ax2.set_xticklabels([])

# plot m
for key in m_results.keys():
    ax0.loglog(np.array(times) / day2sec, m_results[key], label=fr"$M$={float(key) / m_sun:.2f} $M_{{\odot}}$")
    ax0.set_ylabel(r"Luminosity (ergs s$^{-1}$")
    ax0.set_xlabel("time (days)")
    ax0.set_title(r"Varying Ejecta Mass, $M$")
    ax0.legend()

# plot v_sc
for key in v_sc_results.keys():
    ax1.loglog(np.array(times) / day2sec, v_sc_results[key], label=fr"$v_{{sc}}$={float(key):.2E}")
    # ax[1].set_ylabel(r"Luminosity (ergs s$^{-1}$")
    ax1.set_xlabel("time (days)")
    ax1.set_title(r"Varying Velocity Scale, $V_{sc}$")
    ax1.legend()

# plot m_ni
for key in m_ni_results.keys():
    ax2.loglog(np.array(times) / day2sec, m_ni_results[key],
               label=fr"$M_{{Ni}}$={float(key) / m_sun:.2f} $M_{{\odot}}$")
    # ax[2].set_ylabel(r"Luminosity (ergs s$^{-1}$")
    ax2.set_xlabel("time (days)")
    ax2.set_title(r"Varying Nickel Mass, $M_{Ni}$")
    ax2.legend()

plt.tight_layout()
for axis in [ax0.xaxis, ax1.xaxis, ax2.xaxis, ]:
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    axis.set_major_formatter(formatter)
# for axis in [ax0.xaxis, ax0.yaxis, ax1.xaxis, ax1.yaxis, ax2.xaxis, ax2.yaxis,]:

plt.savefig("Parameter_behaviour.png")
plt.show()
