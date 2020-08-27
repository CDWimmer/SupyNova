import tkinter as tk
from tkinter import ttk
from tkinter import N, W, E, S
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from Fitter import _lc_maker
from nagy.nagy import kappa_nagy

day2sec = 86400  # s          | seconds in a day
m_sun = 1.989e33  # g         | solar mass in grams

# todo: some way of loading these automatically or something maybe ?
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

root = tk.Tk()
root.title("SN Data")
root.resizable(True, True)
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

mainframe = ttk.Frame(root, padding="5 5 5 5")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

picker_frame = ttk.Frame(mainframe, padding="5 5 5 5")
picker_frame.grid(column=0, row=1, sticky=(N, W, E, S))

graph_frame = ttk.Frame(mainframe, padding="5 5 5 5")
graph_frame.grid(column=0, row=0, sticky=(N, W, E, S))

fitter_frame = ttk.Frame(mainframe, padding="5 5 5 5")
fitter_frame.grid(column=1, row=0, sticky=(N, W, E, S))

m_value = tk.StringVar()
m_value.set("1.5")
tk.Label(fitter_frame, font=("Helvetica", 16), text="Ejecta Mass (M_sun):").grid(row=0, column=0)
m_input = tk.Entry(fitter_frame, textvariable=m_value)
m_input.grid(row=0, column=1)

mni_value = tk.StringVar()
mni_value.set("0.1")
tk.Label(fitter_frame, font=("Helvetica", 16), text="Nickel Mass (M_sun):").grid(row=1, column=0)
mni_input = tk.Entry(fitter_frame, textvariable=mni_value)
mni_input.grid(row=1, column=1)

vsc_value = tk.StringVar()
vsc_value.set("1e9")
tk.Label(fitter_frame, font=("Helvetica", 16), text="Velocity Scale (cm/s):").grid(row=2, column=0)
vsc_input = tk.Entry(fitter_frame, textvariable=vsc_value)
vsc_input.grid(row=2, column=1)

shift_value = tk.StringVar()
shift_value.set("1.5")
tk.Label(fitter_frame, font=("Helvetica", 16), text="Shift (days):").grid(row=3, column=0)
vsc_input = tk.Entry(fitter_frame, textvariable=shift_value)
vsc_input.grid(row=3, column=1)

data_skip = tk.StringVar()
data_skip.set(2)



chosen = tk.StringVar(root)
chosen.set(files[0])

# graph stuff
fig = plt.Figure(figsize=(8, 6), dpi=100)
ax = fig.add_subplot(111)
plot = FigureCanvasTkAgg(fig, graph_frame)
plot.get_tk_widget().grid(column=0, row=0, sticky=(N, W, E, S))
phase, log_lum, uncert = np.loadtxt(chosen.get(), unpack=True)
time = np.array([(t + abs(phase[0])) * day2sec for t in phase])[int(data_skip.get()):]
log_lum = log_lum[int(data_skip.get()):]
graph, = ax.plot(time/day2sec, log_lum, 'bo', label="Data")
graph2, = ax.plot([1, 2, 3, 4], [1, 2, 3, 4], 'k--', label="Calculated LC")
ax.legend()

# set axis limits to be a bit bigger than this initial plot, because idk what to set.
ax.set_xlim(-time[1]*2/day2sec, time[-1] * 1.2/day2sec)  # these are very far from good solutions
ax.set_ylim(log_lum[0] * 0.9, log_lum[-1] * 1.1)

ax.set_xlabel("time (days)")
ax.set_ylabel(r"log$_{10}$ergs s$^{-1}$")
fig.suptitle(chosen.get())


tk.Label(picker_frame, font=("Helvetica", 17), text="Show data:").grid(row=1, column=0)
tk.Label(picker_frame, font=("Helvetica", 15), text="data skip:").grid(row=1, column=0, sticky=E)
shift_input = tk.Entry(picker_frame, textvariable=data_skip)
shift_input.grid(row=1, column=1)
opt = tk.OptionMenu(picker_frame, chosen, *files)
opt.config(width=50, font=('Helvetica', 12))
opt.grid(column=0, row=2, sticky=(W, E))


alert = tk.Label(fitter_frame, font=("Helvetica", 10), text="Ready.")
alert.grid(row=5, column=0, columnspan=2)


def do_plot(*args):
    times = np.linspace(time[0], time[-1], 300)
    try:
        res = _lc_maker(times, m_sun * float(m_value.get()), float(vsc_value.get()), m_sun * float(mni_value.get()),
                        float(shift_value.get()) * day2sec, kappa_func=kappa_nagy)
        print(res)
        if res.count(np.inf) > 0 or res.count(np.nan) > 0:
            alert['text'] = "Result contains some invalid values!\nSee list printed to console."
        else:
            alert['text'] = 'Ready.'
    except Warning as w:
        print(str(w))
        alert['text'] = f"Warning:\n{str(w)}"
    except Exception as e:
        print(str(e))
        alert['text'] = f"Error:\n{str(e)}"
    else:
        graph2.set_ydata(res)
        graph2.set_xdata(times/day2sec)
        #ax.set_xlim(left=-times[1]/day2sec)
        fig.canvas.draw()
        fig.canvas.flush_events()


def callback(*args):
    phase, log_lum, uncert = np.loadtxt(chosen.get(), unpack=True)
    time = np.array([(t + abs(phase[0])) * day2sec for t in phase])[int(data_skip.get()):]
    log_lum = log_lum[int(data_skip.get()):]
    graph.set_ydata(log_lum)
    graph.set_xdata(time/day2sec)
    ax.set_xlim(left=-time[1]*2/day2sec, right=time[-1]*1.2/day2sec)  # this is very far from a good solution
    ax.set_ylim(log_lum[0] * 0.9, log_lum[-1] * 1.1)  # this is very far from a good solution
    fig.suptitle(chosen.get())
    fig.canvas.draw()
    fig.canvas.flush_events()


root.bind('<Return>', do_plot)

generate = tk.Button(fitter_frame, text="Go", height="2", width="10", command=do_plot)
generate.grid(row=4)

chosen.trace("w", callback)

root.mainloop()
