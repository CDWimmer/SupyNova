import Fitter
from nagy.nagy import kappa_nagy
from arnett_small_R import kappa_const
from random import uniform
m_sun = 1.989e33  # g         | solar mass in grams

files = [
    "data/1998bw_bolom.txt",  # Type Ic
    "data/2009jf_bolom.txt",  # Type Ib
    "data/2011dh_bolom.txt",  # Type II
    "data/2011fe_bolom.txt",  # Type Ia
    #"data/2002ap_nature.txt",
    #"data/2006aj_nature.txt",
    "data/1994I_bolom.txt",  # Type Ic
    "data/1994D_bolom.txt",   # Type Ia
    "data/1991T_bolom.txt",  # Type Ia
    "data/1991T-mod_bolom.txt",  # Type Ia - removed very late time
]

for file in files:
    m_init = 2*m_sun
    vsc_init = 0.9e9
    data_skip = 2
    tries = 0
    while True:
        try:
            Fitter.fit(file, data_skip=data_skip, kappa_function=kappa_nagy,
                       m_init=m_init, vsc_init=vsc_init, m_max=40*m_sun)
        except OverflowError:  # something went bad
            if tries > 10:
                print("Giving up...")
                break
            tries += 1
            print("Something went bad, trying again with randomised initial params")
            m_init = uniform(0.1*m_sun, 10*m_sun)
            vsc_init = uniform(0.01e9, 1e9)
            continue
        else:  # successful fit, move on to next file
            break


for file in files:
    m_init = 7*m_sun
    vsc_init = 1e9
    data_skip = 3
    tries = 0
    while True:
        try:
            Fitter.fit(file, data_skip=data_skip, kappa_function=kappa_const,
                       m_init=m_init, vsc_init=vsc_init, m_max=40*m_sun)
        except OverflowError:  # something went bad
            if tries > 10:
                print("Giving up...")
                break
            tries += 1
            print("Something went bad, trying again with randomised initial params")
            m_init = uniform(0.1*m_sun, 10*m_sun)
            vsc_init = uniform(0.01e9, 1e9)
            continue
        else:  # successful fit, move on to next file
            break
