import tkinter as tk
from tkinter import ttk
from tkinter import N, W, E, S
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np

# todo: some way of loading these automatically or something maybe ?
files = [
    "data/1998bw_bolom.txt",
    "data/2009jf_bolom.txt",
    "data/2011dh_bolom.txt",
    "data/2011fe_bolom.txt",
    "data/2002ap_nature.txt",
    "data/2006aj_nature.txt",
    "data/1994l_bolom.txt",  # Type Ic
    "data/1994D_bolom.txt",   # Type Ia
    "data/1991T_bolom.txt",  # Type Ia
]

root = tk.Tk()
root.title("SN Data")
root.resizable(False, False)
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

mainframe = ttk.Frame(root, padding="5 5 5 5")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))


picker_frame = ttk.Frame(mainframe, padding="5 5 5 5")
picker_frame.grid(column=0, row=0, sticky=(N, W, E, S))

graph_frame = ttk.Frame(mainframe, padding="5 5 5 5")
graph_frame.grid(column=0, row=1, sticky=(N, W, E, S))

chosen = tk.StringVar(root)
chosen.set(files[0])

# graph stuff

fig = plt.Figure(figsize=(6,5), dpi=100)
ax = fig.add_subplot(111)
plot = FigureCanvasTkAgg(fig, graph_frame)
plot.get_tk_widget().grid(column=0, row=0, sticky=(N, W, E, S))
phase, log_lum, uncert = np.loadtxt(chosen.get(), unpack=True)
graph, = ax.plot(phase, log_lum)
ax.set_xlabel("phase (days)")
ax.set_ylabel(r"log$_{10}$ergs s$^{-1}$")
fig.suptitle(chosen.get())


opt = tk.OptionMenu(root, chosen, *files)
opt.config(width=50, font=('Helvetica', 12))
opt.grid(column=0, row=1, sticky=(W, E))


def callback(*args):
    phase, log_lum, uncert = np.loadtxt(chosen.get(), unpack=True)
    graph.set_ydata(log_lum)
    graph.set_xdata(phase)
    fig.suptitle(chosen.get())
    fig.canvas.draw()
    fig.canvas.flush_events()





chosen.trace("w", callback)


root.mainloop()

