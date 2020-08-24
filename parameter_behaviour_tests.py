from arnett_small_R import bolo_l
from nagy.nagy import kappa_nagy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams

# test values:
m_default = 6 * m_sun
m_array = [1.5 * m_sun, 3 * m_sun, 4.5 * m_sun, 6 * m_sun, 7.5 * m_sun, 9 * m_sun]
v_sc_default = 1.2e9
v_sc_array = [0.6e9, 1.2e9, 1.8e9, 2.4e9, 3e9, 3.6e9]
m_ni_default = 0.75 * m_sun
m_ni_array = [0.01 * m_sun, 0.1 * m_sun, 0.25 * m_sun, 0.5 * m_sun, 0.75 * m_sun, 1 * m_sun]
m_ni_array = np.linspace(0.01 * m_sun, m_sun, 6)

skip_days = 1  # Days to skip at start. (shock breakout makes it inacurate)
end_day = 300  # day on which to end simulation
times = np.linspace(skip_days * day2sec, end_day * day2sec, 1000)

m_results = dict()
v_sc_results = dict()
m_ni_results = dict()

# Generate data:

for m in m_array:
    res = []
    for t in times:
        res.append(bolo_l(t, m, v_sc_default, m_ni_default, kappa=kappa_nagy))
    m_results[str(m)] = res
    print(f"For M = {m/m_sun} M_sun, peak of L={max(res):.3E} ergs s^-1 at t={times[np.argmax(res)]/day2sec} days")

for v_sc in v_sc_array:
    res = []
    for t in times:
        res.append(bolo_l(t, m_default, v_sc, m_ni_default, kappa=kappa_nagy))
    v_sc_results[str(v_sc)] = res
    print(f"For v_sc = {v_sc/1e9}e9, peak of L={max(res):.3E} ergs s^-1 at t={times[np.argmax(res)]/day2sec} days")

for m_ni in m_ni_array:
    res = []
    for t in times:
        res.append(bolo_l(t, m_default, v_sc_default, m_ni, kappa=kappa_nagy))
    m_ni_results[str(m_ni)] = res
    print(f"For M_ni = {m_ni/m_sun} M_sun, peak of L={max(res):.3E} ergs s^-1 at t={times[np.argmax(res)]/day2sec} days")


# prepare figures/axes

fig, (ax0, ax1, ax2) = plt.subplots(ncols=3, sharey="all")
fig.set_size_inches(11, 6.4)
ax1.set_xticklabels([])
ax2.set_xticklabels([])

# plot m
for key in m_results.keys():
    ax0.loglog(np.array(times) / day2sec, m_results[key], label=fr"$M$={float(key) / m_sun:.2f} $M_{{\odot}}$")
    ax0.set_ylabel(r"Luminosity (ergs s$^{-1}$)")
    ax0.set_xlabel("time (days)")
    ax0.set_title(r"Varying Ejecta Mass, $M$")
    ax0.legend()

# plot v_sc
for key in v_sc_results.keys():
    ax1.loglog(np.array(times) / day2sec, v_sc_results[key], label=fr"$v_{{sc}}$={float(key)/1e9:.1f}$\times 10^9$")
    # ax[1].set_ylabel(r"Luminosity (ergs s$^{-1}$)")
    ax1.set_xlabel("time (days)")
    ax1.set_title(r"Varying Velocity Scale, $V_{sc}$")
    ax1.legend()

# plot m_ni
for key in m_ni_results.keys():
    ax2.loglog(np.array(times) / day2sec, m_ni_results[key],
               label=fr"$M_{{Ni}}$={float(key) / m_sun:.2f} $M_{{\odot}}$")
    # ax[2].set_ylabel(r"Luminosity (ergs s$^{-1}$)")
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




