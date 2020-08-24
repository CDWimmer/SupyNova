from arnett_small_R import l_ni, bolo_l, l_co
from nagy.nagy import kappa_nagy
from matplotlib import pyplot as plt
import numpy as np

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams

m_array = np.linspace(0.5 * m_sun, 15 * m_sun, 100)
v_sc_array = np.linspace(0.6e8, 10.8e9, 100)
t = 20

mass2lambdani = []
for m in m_array:
    mass2lambdani.append(np.log10(l_ni(t, m, 1.2e9, kappa_nagy)[1]))

mass2lambdaco = []
for m in m_array:
    mass2lambdaco.append(np.log10(l_co(t, m, 1.2e9, kappa_nagy)))

mass2lum = []
for m in m_array:
    mass2lum.append(np.log10(bolo_l(t, m, 1.2e9, 0.75 * m_sun, kappa_nagy)))


vsc2lambdani = []
for v_sc in v_sc_array:
    vsc2lambdani.append(np.log10(l_ni(t, 6 * m_sun, v_sc, kappa_nagy)[1]))

vsc2lambdaco = []
for v_sc in v_sc_array:
    vsc2lambdaco.append(np.log10(l_co(t, 6 * m_sun, v_sc, kappa_nagy)))


vsc2lum = []
for v_sc in v_sc_array:
    vsc2lum.append(np.log10(bolo_l(t, 6 * m_sun, v_sc, 0.75 * m_sun, kappa_nagy)))


fig, ((ax0, ax1, ax2), (ax3, ax4, ax5)) = plt.subplots(ncols=3, nrows=2)
fig.set_size_inches(8.5, 5.2)
ax0.semilogx(m_array / m_sun, mass2lambdani, 'k-')
ax0.set_ylabel(r"$\log_{10}\ \Lambda_{Ni}\left(x, y\right)$")
ax0.set_xlabel(r"$M\ (M_\odot)$")

ax1.semilogx(m_array / m_sun, mass2lambdaco, 'k-')
ax1.set_ylabel(r"$\log_{10}\ \Lambda_{Co}\left(x, y\right)$")
ax1.set_xlabel(r"$M\ (M_\odot)$")

ax2.semilogx(m_array / m_sun, mass2lum, 'k-')
ax2.set_ylabel(r"$\log_{10}\ L_{{bolo}}(t)$")
ax2.set_xlabel(r"$M\ (M_\odot)$")

ax3.semilogx(v_sc_array, vsc2lambdani, 'k-')
ax3.set_ylabel(r"$\log_{10}\ \Lambda\left(x, y\right)$")
ax3.set_xlabel(r"$v_{sc}$")

ax4.semilogx(v_sc_array, vsc2lambdaco, 'k-')
ax4.set_ylabel(r"$\log_{10}\ \Lambda_{Co}\left(x, y\right)$")
ax4.set_xlabel(r"$v_{sc}$")

ax5.semilogx(v_sc_array, vsc2lum, 'k-')
ax5.set_ylabel(r"$\log_{10}\ L_{{bolo}}(t)$")
ax5.set_xlabel(r"$v_{sc}$")


fig.tight_layout()
plt.savefig("Extra_params.png")
plt.show()
