from arnett_small_R import bolo_l, kappa_const
from nagy.nagy import kappa_nagy
from lmfit import Model, Parameters
from collections import Iterable
from numpy import loadtxt, log10, array, sqrt, linspace
from matplotlib import pyplot as plt
from matplotlib.ticker import ScalarFormatter

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams


def fit(data: str, data_skip: int, kappa_function: callable,
        m_min=0.1*m_sun, m_max=30*m_sun, m_init=7*m_sun,
        vsc_min=0.01e9, vsc_max=2e9, vsc_init=1e9,
        mni_min=0.001*m_sun, mni_max=1.5 * m_sun, mni_init=0.1 * m_sun):

    # ============== EDIT THESE ============== #
    # data = "data/2009jf_bolom.txt"
    # data_skip = 4
    # below should be an instance of a function that takes 1 arg: time, t.
    # kappa_function = kappa_const
    # ======================================== #

    def fit_me(times: Iterable, m, v_sc, m_ni, shift):
        """
        Return bolometric luminosity values for a given set of times.
        Uses the kappa_nagy function for kappa calculation.
        """
        time_shifted = [t + shift for t in times]
        results = []
        for t in time_shifted:
            results.append(log10(bolo_l(t, m, v_sc, m_ni, kappa=kappa_function)))
        return results

    # load data file
    phase, log_lum, uncert = loadtxt(data, unpack=True)

    time = array([(t + abs(phase[0])) * day2sec for t in phase])[data_skip:]
    log_lum = log_lum[data_skip:]


    lc_model = Model(fit_me)
    print('parameter names: {}'.format(lc_model.param_names))
    print('independent variables: {}'.format(lc_model.independent_vars))

    params = Parameters()
    # ============== EDIT THESE? ============= #
    params.add(name="m", value=m_init, min=m_min, max=m_max)
    params.add(name="v_sc", value=vsc_init, min=vsc_min, max=vsc_max)
    params.add(name="mni_minus_m", value=1, vary=True, min=1e-12)
    # params.add(name="m_ni", value=mni_init, min=mni_min, max=mni_max)
    params.add(name="m_ni", expr='0.2*m + mni_minus_m')  # max value of m_ni is 0.2*m
    params.add(name="shift", value=500, min=1, max=15 * day2sec)
    # ======================================== #

    # Do a fit:
    print(f"\033[94mFitting {data}....\033[0m")
    try:
        result = lc_model.fit(log_lum, params, times=time)
    except Exception as e:
        raise OverflowError(f"Something went wrong during fitting, perhaps parameters went mad, try different initial "
                            f"values?\n\033[91m{str(e)}\033[0m")
        # print("Showing plot of the data and initial params")
        # print(type(time), type(log_lum))
        # plt.scatter(time/day2sec, log_lum, 'bo', label="data")
        # plt.plot(time/day2sec, fit_me(time, params["m"].value, params["v_sc"].value, params["m_ni"].value, 500), 'k--',
        #          label="init params")
        # plt.legend()
        # plt.show()
    else:  # it didn't crash
        print("Done!")
        # print(result.fit_report())  # dumps a load of fitting data

        # generate a smoother plot than result.best_fit provides:
        # extract fit results
        m_result = result.params['m'].value
        v_sc_result = result.params['v_sc'].value
        m_ni_result = result.params['m_ni'].value
        shift_result = result.params["shift"].value
        chi2_result = result.chisqr
        red_chi2_result = result.redchi

        smooth_times = linspace(time[0], time[-1], len(time) * 3)  # 3 times as many time points, and evenly spaced
        init_result = fit_me(  # smooth curve of initial guess
            smooth_times,
            result.init_values["m"],
            result.init_values["v_sc"],
            result.init_values["m_ni"],
            result.init_values["shift"])
        best_result = fit_me(  # smooth curve of best guess
            smooth_times,
            m_result,
            v_sc_result,
            m_ni_result,
            shift_result)

        fig, ax = plt.subplots(1, 1)
        ax.semilogx(time / day2sec, log_lum, 'bo', label=f"Data for {data.split('/')[1].split('.')[0].split('_')[0]}")
        ax.plot(smooth_times / day2sec, init_result, 'k--', label="Initial fit")
        ax.plot(smooth_times / day2sec, best_result, 'r-', label="Best fit")
        ax.set_xlabel("time (days)")
        ax.set_ylabel(r"log$_{10}$ luminosity")
        ax.legend()
        formatter = ScalarFormatter()
        formatter.set_scientific(False)
        ax.xaxis.set_major_formatter(formatter)
        result_string = (f"kappa function: {kappa_function.__name__}\n"
                         f"Data skip: {data_skip}\n"
                         f"M = {m_result / m_sun:.5f} solar masses\n"
                         f"v_sc = {v_sc_result / 1e9:.5f} *10^9 cm/s\n"
                         f"M_Ni = {m_ni_result / m_sun:.5f} solar masses\n"
                         f"shift = {shift_result/day2sec:.5f} days\n"
                         f"chi squared = {chi2_result}\n"
                         f"reduced chi = {red_chi2_result}\n"
                         f"sqrt(2.53M_Ni/M) = {sqrt(2.53 * m_ni_result / m_result)}")
        print(result_string)

        print("Writing results to disk...")
        try:
            with open(f"fits/{data.split('/')[1].split('.')[0]}_{kappa_function.__name__}.txt", 'w') as f:
                f.write(result_string)
            fig.savefig(f"fits\\{data.split('/')[1].split('.')[0]}_{kappa_function.__name__}.png")
        except Exception as e:
            print("Something went wrong when saving!", str(e))
        else:
            print("Done!")  # Showing plot.")
        # plt.show()

# todo: make this usable from command line?
