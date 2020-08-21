from arnett_small_R import l_ni, l_co, bolo_l, kappa_const
import numpy as np
from nagy import nagy
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

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


skip_days = 1  # Days to skip at start. (shock breakout makes it inacurate)
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
    lum_res.append(bolo_l(t, m=m, v_sc=v_sc, m_ni=m_ni))

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
    lum_res_var.append(bolo_l(t, m=m, v_sc=v_sc, m_ni=m_ni, kappa=nagy.kappa_nagy))

# plotting
print(f"Using test values of m=1.45*m_sun g, v_sc=1.2e9 cm/s ")
fig, ax = plt.subplots(1, 2)
fig.set_size_inches(10, 5.4)

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

for axis in [ax[0].xaxis, ax[1].xaxis,]:
    formatter = ScalarFormatter()
    formatter.set_scientific(False)
    axis.set_major_formatter(formatter)

plt.tight_layout()
plt.savefig("demo.png")
plt.show()
