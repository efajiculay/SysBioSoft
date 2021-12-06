import sys
import os

from datetime import datetime
import tkinter as gui
from tkinter import ttk
from tkinter import filedialog
from pathlib import Path
from PIL import Image as Image2, ImageTk
import numpy as np
import time
import threading
from queue import Queue
from math import ceil as Myceil
from sys import platform
import webbrowser

if platform == "win32":
    from subprocess import Popen, CREATE_NEW_CONSOLE
elif platform == "darwin":
    from applescript import tell as my_tell_us
    from subprocess import check_output, call as my_call_us
elif platform == "linux":
    pass
else:
    from subprocess import Popen

try:
    import tempfile
    Temporary_folder = str(tempfile.gettempdir()).replace(
        "\\", "/")+"/BioSANS_temporary_folder"
except:
    Temporary_folder = "BioSANS_temporary_folder"

from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global as proc_global
from BioSANS2020.gui_functs.prepare_canvas import *
from BioSANS2020.prepcodes.process import *
from BioSANS2020.analysis.plotting.plot_traj import *
from BioSANS2020.analysis.numeric.transform_data import *
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.model.ode_parse import odeExtract

import BioSANS2020.model.topology_view as topology_view
from BioSANS2020.model.new_file import *

globals2.init()
if __name__ == '__main__':
    proc_global.init()

top = gui.Tk()
top.title("BioSANS 1.0")
top.geometry("1005x550")
#top.resizable(False, False)

header = gui.Label(top, text="BioSANS")
header.configure(
    bg="green",
    fg="white",
    height=1,
    #width = 1005,
    font="Helvetica 18 bold italic"
)
header.pack(fill="x")

frame = gui.Frame(top)
frame.configure(
    bg="light cyan",
    borderwidth=2,
    height=500,
    width=1005
)
frame.pack(fill="both", expand=True)

footer = gui.Label(
    top, text="Biological Symbolic and Numeric Simulation Algorithms")
footer.configure(
    bg="green",
    fg="white",
    #width = 1005,
    font="Helvetica 10 bold italic",
    anchor='w'
)
footer.pack(fill="x")

file_name = {"topology": Temporary_folder, "current_folder": None}


def load_data(items):
    global file_name
    file = filedialog.askopenfilename(title="Select file")
    file_name["topology"] = file
    file_name["current_folder"] = file
    globals2.toConvert = file
    if os.path.isfile(file):
        file_name['last_open'] = topology_view.view_topo(file, items)


def show_file_dir(path):
    global platform
    if platform == "win32":
        os.startfile(os.path.dirname(path))
    elif platform == "darwin":
        my_call_us('open', os.path.dirname(path))
    elif platform == "linux":
        my_call_us('xdg-open', os.path.dirname(path))
    else:
        webbrowser.open(os.path.dirname(path))


def create_file(items, Ftype):
    global file_name, Temporary_folder
    try:
        os.mkdir(Temporary_folder, 0o777)
    except:
        for item in Path(Temporary_folder).iterdir():
            if item.is_dir():
                pass
            else:
                item.unlink()

    file_name['last_open'] = new_file(items)
    file_name["topology"] = Temporary_folder+"/temp.txt"
    if Ftype == 1:
        file_name['last_open'].insert(INSERT, "FUNCTION_DEFINITIONS:\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(
            INSERT, "#REACTIONS, Volume = 1, tend = 100, steps = 100, FileUnit = molar\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "@CONCENTRATION\n")
    elif Ftype == 2:
        file_name['last_open'].insert(INSERT, "ODE_DECLARATIONS:\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "INI_CONCENTRATIONS:\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "\n")
        file_name['last_open'].insert(INSERT, "RATE_CONSTANTS:\n")


def extractODE(items):
    global file_name
    file_name['last_open'] = odeExtract.odedxdt_to_topo(
        file_name["topology"], items)
    file_name["topology"] = file_name["topology"]+".top"


def sbml_to_topo2(tocon, items):
    file_name["topology"] = sbml_to_topo(tocon)
    file_name['last_open'] = topology_view.view_topo(
        file_name["topology"], items)


def save_file():
    global file_name
    file = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    file_name["topology"] = file.name
    file_name["current_folder"] = file.name
    if file is None:
        return
    file.write(file_name['last_open'].get("0.0", END))
    file.close()


def runpy_file():
    global file_name, PIPE
    with open(file_name["topology"], "w") as ffvar:
        ffvar.write(file_name['last_open'].get("0.0", END))
        ffvar.write("\ninput('Press enter to exit:')")
    if platform == "win32":
        Popen([sys.executable, file_name["topology"]],
              creationflags=CREATE_NEW_CONSOLE)
    elif platform == "darwin":
        my_tell_us.app("Terminal", 'do script "' +
                       str(sys.executable)+" "+file_name["topology"]+'"')
    elif platform == "linux":
        os.system("gnome-terminal -x python3 "+file_name["topology"])
    else:
        Popen(str(sys.executable)+" "+file_name["topology"], shell=True)


def run_SSL():
    if platform == "win32":
        Popen([sys.executable, "-m", "BioSANS2020.BioSSL"],
              creationflags=CREATE_NEW_CONSOLE)
    elif platform == "darwin":
        A = check_output(['pip3', 'show', 'BioSANS2020-efajiculay'])
        A = str(A).split("\\n")
        install_dir = ""
        for row in A:
            line = row.split(":")
            if line[0].strip() == "Location":
                install_dir = ("".join(line[1:]).strip("\\r\\n").replace(
                    "c\\", "c:/").replace("\\", "/").replace("//", "/"))
        if install_dir != "":
            install_dir = str(install_dir)
            my_tell_us.app("Terminal", 'do script "'+str(sys.executable) +
                           " "+install_dir+"/BioSANS2020/BioSSL.py"+'"')
    elif platform == "linux":
        os.system("gnome-terminal -x python3 " +
                  os.path.join(os.getcwd(), "BioSSL.py"))
    else:
        Popen(str(sys.executable)+" " +
              os.path.join(os.getcwd(), "BioSSL.py"), shell=True)


def load_data2(plot=False):
    t_o = time.time()
    global file_name, current_data
    file = filedialog.askopenfilename(title="Select file")
    file_name["trajectory"] = file
    file_name["current_folder"] = file
    with open(file_name["trajectory"], "r") as f:
        data = []
        ddvar = []
        row1 = str(f.readline()).strip()
        slabels = row1.split("\t")[1:]
        for row in f:
            cols = [float(x) for x in row.split("\t")]
            if cols[0] == 0.0 and len(ddvar) > 0:
                data.append(np.array(ddvar))
                ddvar = []
            ddvar.append(cols)
        data.append(np.array(ddvar))
    if plot:
        plot_traj(data, slabels, items, globals2.plotted, mix_plot=True,
                  logx=False, logy=False, normalize=False)
    print(data[0],"\n\n")
    print(data[1],"\n\n")
    print(data[2],"\n")
    current_data = (data, slabels)
    gui.messagebox.showinfo("showinfo", "Trajectory loaded succesfully")
    # print(time.time()-t_o)


def tload_data2(plot=False):
    if __name__ == '__main__':
        t = threading.Thread(target=load_data2, args=(plot,), daemon=False)
        t.start()


def delete_this(frame, canvas):
    canvas.delete(frame)
    canvas_update_widgets(None, canvas)


def canvas_update_widgets(e, canvas):
    #R = canvas._root().winfo_height()
    #root = canvas._root()
    #width = root.winfo_screenwidth()
    #height = root.winfo_screenheight()
    # print(canvas.winfo_children())
    H = canvas.winfo_height()
    W = canvas.winfo_width()
    objCan = canvas.find_all()
    objLen = len(objCan)
    ind = 0
    for x in objCan:
        canvas.itemconfig(x, height=H, width=W-5)
        x1, y1, x2, y2 = canvas.bbox(x)
        canvas.move(x, 0, (objLen-ind)*H-y1+3)
        ind = ind + 1
    canvas.configure(scrollregion=canvas.bbox("all"))
    return "break"


def load_image(wdata=False):
    t_o = time.time()
    global items, current_data
    canvas, scroll_x, scroll_y = items
    file = filedialog.askopenfilename(title="Select file")
    file_name["current_folder"] = file
    load = Image2.open(file)
    render = ImageTk.PhotoImage(load)
    img = gui.Label(canvas, image=render)
    img.image = render

    fframe = canvas.create_window(
        0, 426*globals2.plot_i, anchor='nw', window=img)
    B = Button(img, text=" X ", fg='red', highlightcolor='blue', bg='white',
               height=1, relief='raised', command=lambda: delete_this(fframe, canvas))
    B.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")

    canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
    canvas.configure(scrollregion=canvas.bbox("all"))
    canvas.bind("<Configure>", lambda e: canvas_update_widgets(e, canvas))
    globals2.plot_i = globals2.plot_i+1

    if wdata:
        file = str(file).replace("jpg", "dat").replace("png", "dat")
        with open(file, "r") as f:
            data = []
            ddvar = []
            row1 = str(f.readline()).strip()
            slabels = row1.split("\t")[1:]
            for row in f:
                cols = [float(x) for x in row.split("\t")]
                if cols[0] == 0.0 and len(ddvar) > 0:
                    data.append(np.array(ddvar))
                    ddvar = []
                ddvar.append(cols)
            data.append(np.array(ddvar))
        current_data = (data, slabels)
    print(time.time()-t_o)
    canvas_update_widgets(None, canvas)


def eval2(x):
    try:
        return eval(x)
    except:
        try:
            par = str(x).lower().capitalize()
            return eval(par)
        except:
            return str(x)


current_data = None


def dict_trans(x1):
    x2 = x1.split(",")
    x3 = {}
    try:
        for x in x2:
            r = x.split("=")
            x3[r[0].strip()] = float(r[1])
    except:
        pass
    return x3


def convert(x, con):
    try:
        return con(x)
    except:
        return x


def range_trans(x1):
    x3 = []
    x2 = x1.split(",")
    if len(x2) > 1:
        x3 = x2[1:]
        x3.append(x2[0])
    return x3


def range_prep(x1):
    x3 = []
    x2 = x1.split(",")
    if len(x2) > 1:
        x2 = [convert(x2[x], float) for x in range(len(x2))]
        cc = x2[0].lower().split("f")
        r2 = 0
        if len(cc) < 2:
            cc = x2[0].lower().split("b")
            r2 = 1
        r1 = int(cc[1])-1
        cc = ",".join([str(x) for x in x2[1:]])
        x3 = list(eval(cc))
        return [[r1, r2], x3]
    return x1


def mrun_propagation(par, E, defs2):
    global items, current_data, file_name
    try:
        del file_name["trajectory"]
    except:
        pass

    defs = []
    for i in range(len(E)):
        try:
            val = eval2(E[i].get())
            defs.append(val)
        except:
            val = eval2(defs2[i].get())
            defs.append(val)
    #defs[15] = dict_trans(E[15].get())
    defs[16] = range_trans(E[16].get())
    defs[17] = range_prep(E[17].get())
    if not defs[9] in ["Analyt", "SAnalyt", "Analyt-ftx", "SAnalyt-ftk", "k_est1", "k_est2", "k_est3", "k_est4", "k_est5", "k_est6", "k_est7", "k_est8", "k_est9", "k_est10", "k_est11", "NetLoc1", "NetLoc2"]:
        with open(defs[13]+"_"+defs[9]+"_params.dat", "w") as f:
            f.write("\n".join([str(x) for x in defs]))
    if defs[9] in ["k_est1", "k_est2", "k_est3", "k_est4", "k_est6", "k_est7", "k_est8", "k_est9", "k_est10", "k_est11", "LNA2", "LNA3", "Analyt", "Analyt-ftx", "SAnalyt", "SAnalyt-ftk", "Analyt2", "topoTosbml", "topoTosbml2", "topoTosbml3", "LNA-vs", "LNA-ks", "LNA-xo", "NetLoc1", "NetLoc2"]:
        par.destroy()
    defs.append(items)
    current_data = tprocess(defs)


super_thread_run = None


def tprocess(defs):
    global super_thread_run
    if defs[9] == "k_est5":
        process(*defs)
    else:
        if __name__ == '__main__':
            super_thread_run = 1
            t = threading.Thread(target=lambda: process(*defs), daemon=False)
            t.start()


def analysis_case(anaCase, items):
    global current_data, file_name
    if "trajectory" not in file_name:
        gui.messagebox.showinfo("Trajectory not loaded yet",
                                "Please load the trajectory. BioSANS save it into a file during your last run.")
        try:
            load_data2(False)
        except:
            try:
                del file_name["trajectory"]
            except:
                pass
            return None
        data, slabels = current_data

    if anaCase == "cov":
        return calc_covariance(current_data, items)
    elif anaCase == "fanoF":
        return fano_factor(current_data, items)
    elif anaCase == "corr":
        return calc_cross_corr(current_data, items)
    elif anaCase == "pdens1":
        return prob_density_calc(current_data, items)
    elif anaCase == "pdens2":
        return prob_density_calc2(current_data, items)
    elif anaCase == "pdens3":
        return prob_density_calc3(current_data, items)
    elif anaCase == "avetrj":
        return ave_traj_calc(current_data, items)
    elif anaCase == "phaseP":
        return plot_trajD(current_data, items)
    elif anaCase == "plotD":
        return plot_trajD2(current_data, items)


def plot_trajD(current_data, items):
    global super_thread_run
    if super_thread_run == 1:
        gui.messagebox.showinfo(
            "Warning: Unsafe thread", "If you want to have a 3D Phase portrait, restart BioSANS and load trajectory. For 2D phase portrait, just continue")
    try:
        data, slabels = current_data
        par = gui.Toplevel()
        par.resizable(False, False)
        par.wm_title("Plot settings")

        L = [
            gui.Label(par, text="choose x-axis", fg="blue"),
            gui.Label(par, text="choose y-axis", fg="blue"),
            gui.Label(par, text="choose z-axis", fg="blue"),
            gui.Label(par, text="choose step-range", fg="blue")
        ]
        [L[i].grid(row=i, column=0, sticky=gui.W, pady=2) for i in range(4)]
        optsel = [[par, gui.StringVar(), "None", "time"] +
                  slabels for i in range(3)]
        optselVar = [optsel[i][1] for i in range(3)]
        [optselVar[i].set("None") for i in range(3)]
        E = [gui.OptionMenu(*optsel[i]) for i in range(3)]
        [E[i].config(width=14) for i in range(3)]
        E.append(gui.Entry(par, bd=5))
        [E[i].grid(row=i, column=1, sticky=gui.W, pady=2) for i in range(4)]
        E[-1].insert(gui.END, "0:-1")

        B1 = ttk.Button(par, text="PLOT", command=lambda: plot_traj2(
            data, slabels, items, globals2.plotted, logx=False, logy=False, normalize=False,
            xlabel=optselVar[0].get(), ylabel=optselVar[1].get(), zlabel=optselVar[2].get(), trange=E[-1].get()))
        B1.grid(row=5, column=0, sticky=gui.W, pady=2)
    except:
        gui.messagebox.showinfo("Trajectory not loaded yet",
                                "Please load the trajectory. BioSANS save it into a file during your last run.")


def getChecked(L1, L, slabels):
    checkSi = []
    for i in range(len(L1)):
        key = L1[i].get()
        if key != "0":
            checkSi.append(slabels.index(key))
    return checkSi


def plot_trajD2(current_data, items):
    try:
        data, slabels = current_data
        pard = gui.Toplevel()
        pard.resizable(True, True)
        pard.wm_title("Plot species")
        pard.maxsize(width=300, height=700)
        canvas, scroll_x, scroll_y = prepare_frame_for_plot(
            pard, width=200, height=200)
        par = gui.Frame(canvas, width=200, height=int(len(slabels)*150))
        par.pack(side="left", fill="both", expand=True)
        pard.configure(
            bg="light blue",
            borderwidth=2,
        )
        canvas.configure(
            bg="light blue",
            borderwidth=2,
        )

        Ls = len(slabels)
        L1 = [gui.StringVar() for i in range(Ls)]
        [L1[i].set("0") for i in range(Ls)]
        L = [gui.Checkbutton(par, text=slabels[i], variable=L1[i],
                             onvalue=slabels[i], offvalue="0") for i in range(Ls)]
        [L[i].grid(row=i, column=0, sticky=gui.W, pady=2)
         for i in range(Myceil(Ls/2))]
        [L[i].grid(row=i-Myceil(Ls/2), column=1, sticky=gui.W, pady=2)
         for i in range(Myceil(Ls/2), Ls)]

        canvas.create_window(0, 0, anchor='n', window=par)
        canvas.configure(yscrollcommand=scroll_y.set,
                         xscrollcommand=scroll_x.set)
        canvas.configure(scrollregion=canvas.bbox("all"))

        B1 = ttk.Button(pard, text="PLOT",
                        command=lambda:
                        plot_traj(data, slabels, items, globals2.plotted, mix_plot=True, logx=False,
                                  logy=False, normalize=False, si_ticked=getChecked(L1, L, slabels))
                        )
        B1.pack(side="bottom", fill="x")
    except:
        gui.messagebox.showinfo("Trajectory not loaded yet",
                                "Please load the trajectory. BioSANS save it into a file during your last run.")


def paramSet(method):
    global file_name, items
    with open(file_name["topology"], "w") as ffvar:
        ffvar.write(file_name['last_open'].get("0.0", END))

    path = Path(file_name["topology"])
    ss = str(file_name["topology"]).split("/")
    ss = ss[-1] if len(ss) > 1 else ""
    name = str(path.parent)+"/"+ss+"_"+datetime.now().strftime("%Y%m%d_%H%M%S")
    par = gui.Toplevel()
    par.resizable(False, False)
    par.wm_title("Parameter setting")
    opts = [
        "File name",
        "Number of iteration :",
        "File Units? :",
        "Volume (L) :",
        "end time (tn) :",
        "tau-scaler",
        "Normalized",
        "logx",
        "logy",
        "method",
        "tsteps",
        "mix_plot",
        "save",
        "out fname",
        "show plot",
        "time label",
        "Cini range",
        "K-range",
        "mult proc",
        "Implicit"
    ]
    defs = [file_name["topology"], 1, gui.StringVar(), 1.0, 100, 1.5, False, False,
            False, method, 1000, True, True, name, True, "time (sec)", "", "", False, False]
    defs[2].set('molecules')

    topfile = open(file_name["topology"], "r")
    for row in topfile:
        if row[0] == "#":
            gg = row.split(",")[1:]
            for x in gg:
                xx = [g.strip() for g in x.split("=")]
                if xx[0] == "Volume":
                    defs[3] = xx[1]
                elif xx[0] == "tend":
                    defs[4] = xx[1]
                elif xx[0] == "FileUnit":
                    defs[2].set(xx[1])
                elif xx[0] == "logx":
                    defs[7] = xx[1]
                elif xx[0] == "Normalized":
                    defs[6] = xx[1]
                elif xx[0] == "steps":
                    defs[10] = xx[1]

    oplen = len(opts)
    L = [gui.Label(par, text=opts[i], fg="blue") for i in range(oplen)]
    [L[i].grid(row=i, column=0, sticky=gui.W, pady=2) for i in range(10)]
    [L[10+i].grid(row=i, column=2, sticky=gui.W, pady=2) for i in range(10)]
    E = [gui.Entry(par, bd=5) if i != 2 else gui.OptionMenu(
        par, defs[2], 'molecules', 'molar', 'moles') for i in range(oplen)]
    [E[i].grid(row=i, column=1, sticky=gui.W, pady=2) for i in range(10)]
    [E[10+i].grid(row=i, column=3, sticky=gui.W, pady=2) for i in range(10)]
    [E[i].insert(gui.END, str(defs[i])) if i !=
     2 else None for i in range(oplen)]
    E[9].configure(state="disable")
    E[2].config(width=14)
    if method == "ODE":
        E[5].configure(state="disable")
    if method == "ODE2":
        E[5].configure(state="disable")
    elif method == "Gillespie_":
        E[5].configure(state="disable")
    elif method == "CLE":
        E[5].delete(0, gui.END)
        E[5].insert(gui.END, str(10))

    B1 = ttk.Button(par, text="RUN",
                    command=lambda: mrun_propagation(par, E, defs))
    B1.grid(row=oplen, column=0, sticky=gui.W, pady=2)
    if method in ["k_est1", "k_est2", "k_est3", "k_est4", "k_est6", "k_est7", "k_est8", "k_est9", "k_est10", "k_est11", "LNA2", "LNA3", "LNA-vs", "LNA-ks", "LNA-xo", "NetLoc1", "NetLoc2", "Analyt", "SAnalyt-ftk", "SAnalyt", "Analyt-ftx", "Analyt2", "topoTosbml", "topoTosbml2", "topoTosbml3"]:
        B1.invoke()


if __name__ == "__main__":

    menubut1 = gui.Menubutton(frame, text=" File/Model ", activebackground="#f2f20d",
                              activeforeground="red", bg="#00cc00", fg="white" if platform.lower() != "darwin" else "green")
    menubut1.menu = gui.Menu(menubut1, tearoff=1)
    menubut1["menu"] = menubut1.menu
    LoadMenu = gui.Menu(frame, tearoff=1)
    LoadMenu.add_command(label="Topology/File", command=lambda: load_data(items),
                         background="white", foreground="Blue")
    LoadMenu.add_command(label="Trajectory file",
                         command=tload_data2, background="white", foreground="Blue")
    LoadMenu.add_command(label="Traj. w/ plot", command=lambda: tload_data2(True),
                         background="white", foreground="Blue")
    LoadMenu.add_command(label="Image of plot", command=load_image,
                         background="white", foreground="Blue")
    LoadMenu.add_command(label="Image w/ data", command=lambda: load_image(True),
                         background="white", foreground="Blue")
    LoadMenu.add_command(label="Current folder", command=lambda: show_file_dir(
        file_name["current_folder"]), background="white", foreground="Blue")
    menubut1.menu.add_cascade(label="Open", menu=LoadMenu)
    NewFMenu = gui.Menu(frame, tearoff=1)
    NewFMenu.add_command(label="Blank File", command=lambda: create_file(
        items, 0), background="white", foreground="Blue")
    NewFMenu.add_command(label="Topo File", command=lambda: create_file(
        items, 1), background="white", foreground="Blue")
    NewFMenu.add_command(label="ODE File", command=lambda: create_file(
        items, 2), background="white", foreground="Blue")
    menubut1.menu.add_cascade(label="New File", menu=NewFMenu)
    menubut1.menu.add_command(label="Save File", command=lambda: save_file())
    menubut1.menu.add_command(
        label="Run File.py", command=lambda: runpy_file())
    menubut1.menu.add_command(label="Run SSL", command=lambda: run_SSL())
    ConvMenu = gui.Menu(frame, tearoff=1)
    ConvMenu.add_command(label="SBML to Topo", command=lambda: sbml_to_topo2(
        globals2.toConvert, items), background="white", foreground="Blue")
    ConvMenu.add_command(label="ODE to Topo", command=lambda: extractODE(
        items), background="white", foreground="Blue")
    TopSbml = gui.Menu(frame, tearoff=1)
    TopSbml.add_command(label="molecules", command=lambda: paramSet(
        "topoTosbml"), background="white", foreground="Blue")
    TopSbml.add_command(label="molar", command=lambda: paramSet(
        "topoTosbml2"), background="white", foreground="Blue")
    TopSbml.add_command(label="no unit", command=lambda: paramSet(
        "topoTosbml3"), background="white", foreground="Blue")
    ConvMenu.add_cascade(label="Topo to SBML", menu=TopSbml)
    menubut1.menu.add_cascade(label="Convert model", menu=ConvMenu)
    ParEsMenu = gui.Menu(frame, tearoff=1)
    ParEsMenu.add_command(label="Nelder-Mead (NM), Macroscopic",
                          command=lambda: paramSet("k_est6"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="Nelder-Mead (NM), Microscopic",
                          command=lambda: paramSet("k_est7"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="Powell, Macroscopic", command=lambda: paramSet(
        "k_est8"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="Powell, Microscopic", command=lambda: paramSet(
        "k_est9"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="L-BFGS-B, Macroscopic", command=lambda: paramSet(
        "k_est10"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="L-BFGS-B, Microscopic", command=lambda: paramSet(
        "k_est11"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="NM-Diff. Evol., Macroscopic",
                          command=lambda: paramSet("k_est3"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="NM-Diff. Evol., Microscopic",
                          command=lambda: paramSet("k_est4"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="Parameter slider/scanner",
                          command=lambda: paramSet("k_est5"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="MCEM, Macroscopic", command=lambda: paramSet(
        "k_est1"), background="white", foreground="Blue")
    ParEsMenu.add_command(label="MCEM, Microscopic", command=lambda: paramSet(
        "k_est2"), background="white", foreground="Blue")
    menubut1.menu.add_cascade(label="Estimate Params", menu=ParEsMenu)
    menubut1.place(x=2, y=5)

    menubut2 = gui.Menubutton(frame, text="Propagation", activebackground="#f2f20d",
                              activeforeground="red", bg="#00cc00", fg="white" if platform.lower() != "darwin" else "green")
    menubut2.menu = gui.Menu(menubut2, tearoff=1)
    menubut2["menu"] = menubut2.menu
    AnalMenu = gui.Menu(frame, tearoff=1)
    AnalMenu.add_command(label="Pure Symbolic :f(t,xo,k)", command=lambda: paramSet(
        "Analyt"), background="white", foreground="Blue")
    AnalMenu.add_command(label="Semi-Symbolic :f(t)", command=lambda: paramSet(
        "SAnalyt"), background="white", foreground="Blue")
    AnalMenu.add_command(label="Semi-Symbolic :f(t,xo)", command=lambda: paramSet(
        "Analyt-ftx"), background="white", foreground="Blue")
    AnalMenu.add_command(label="Semi-Symbolic :f(t,k)", command=lambda: paramSet(
        "SAnalyt-ftk"), background="white", foreground="Blue")
    AnalMenu.add_command(label="For wxmaxima", command=lambda: paramSet(
        "Analyt2"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="Analytical soln.", menu=AnalMenu)

    Rungek4 = gui.Menu(frame, tearoff=1)
    Rungek4.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "rk4-1"), background="white", foreground="Blue")
    Rungek4.add_command(label="Molar(macro)", command=lambda: paramSet(
        "rk4-2"), background="white", foreground="Blue")
    Rungek4.add_command(label="Mole(macro)", command=lambda: paramSet(
        "rk4-3"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="RK4-fix-interval", menu=Rungek4)

    Rungek4a = gui.Menu(frame, tearoff=1)
    Rungek4a.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "rk4-1a"), background="white", foreground="Blue")
    Rungek4a.add_command(label="Molar(macro)", command=lambda: paramSet(
        "rk4-2a"), background="white", foreground="Blue")
    Rungek4a.add_command(label="Mole(macro)", command=lambda: paramSet(
        "rk4-3a"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="RK4-tau-adaptive", menu=Rungek4a)

    EulrTau = gui.Menu(frame, tearoff=1)
    EulrTau.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "Euler-1"), background="white", foreground="Blue")
    EulrTau.add_command(label="Molar(macro)", command=lambda: paramSet(
        "Euler-2"), background="white", foreground="Blue")
    EulrTau.add_command(label="Mole(macro)", command=lambda: paramSet(
        "Euler-3"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="Euler (tau-adaptive-1)", menu=EulrTau)

    EulrTau2 = gui.Menu(frame, tearoff=1)
    EulrTau2.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "Euler2-1"), background="white", foreground="Blue")
    EulrTau2.add_command(label="Molar(macro)", command=lambda: paramSet(
        "Euler2-2"), background="white", foreground="Blue")
    EulrTau2.add_command(label="Mole(macro)", command=lambda: paramSet(
        "Euler2-3"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="Euler (tau-adaptive-2)", menu=EulrTau2)

    ODEMenu = gui.Menu(frame, tearoff=1)
    ODEMenu.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "ODE-1"), background="white", foreground="Blue")
    ODEMenu.add_command(label="Molar(macro)", command=lambda: paramSet(
        "ODE-2"), background="white", foreground="Blue")
    ODEMenu.add_command(label="Mole(macro)", command=lambda: paramSet(
        "ODE-3"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="ODE int", menu=ODEMenu)

    #Itoints = gui.Menu(frame,tearoff = 1 )
    #Itoints.add_command ( label="Molecules",command=lambda: paramSet("Itoint-1"),background="white",foreground="Blue" )
    #Itoints.add_command ( label="Molar",command=lambda: paramSet("Itoint-2"),background="white",foreground="Blue" )
    #Itoints.add_command ( label="Mole",command=lambda: paramSet("Itoint-3"),background="white",foreground="Blue" )
    #menubut2.menu.add_cascade(label="Itoint", menu=Itoints)

    #Stratint = gui.Menu(frame,tearoff = 1 )
    #Stratint.add_command ( label="Molecules",command=lambda: paramSet("Stratint-1"),background="white",foreground="Blue" )
    #Stratint.add_command ( label="Molar",command=lambda: paramSet("Stratint-2"),background="white",foreground="Blue" )
    #Stratint.add_command ( label="Mole",command=lambda: paramSet("Stratint-3"),background="white",foreground="Blue" )
    #menubut2.menu.add_cascade(label="Stratint", menu=Stratint)

    cletauA = gui.Menu(frame, tearoff=1)
    cletauA.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "CLE"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="CLE (tau-adaptive)", menu=cletauA)

    cletauA2 = gui.Menu(frame, tearoff=1)
    cletauA2.add_command(label="Molecules(micro)", command=lambda: paramSet(
        "CLE2"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="CLE (cle-fixIntvl)", menu=cletauA2)

    TauLMenu = gui.Menu(frame, tearoff=1)
    TauLMenu.add_command(label="Tau-leapingV1-micro", command=lambda: paramSet(
        "Tau-leaping"), background="white", foreground="Blue")
    TauLMenu.add_command(label="Tau-leapingV2-micro", command=lambda: paramSet(
        "Tau-leaping2"), background="white", foreground="Blue")
    TauLMenu.add_command(label="Sim-TauLeap-micro", command=lambda: paramSet(
        "Sim-TauLeap"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="Tau-leaping-micro", menu=TauLMenu)

    LNAMenu = gui.Menu(frame, tearoff=1)
    LNAMenu.add_command(label="COV-time-dependent, Macroscopic",
                        command=lambda: paramSet("LNA(t)"), background="white", foreground="Blue")
    LNAMenu.add_command(label="FF-time-dependent, Macroscopic",
                        command=lambda: paramSet("LNA2(t)"), background="white", foreground="Blue")
    LNAMenu.add_command(label="Numeric, values", command=lambda: paramSet(
        "LNA"), background="white", foreground="Blue")
    LNAMenu.add_command(label="Symbolic, Microscopic", command=lambda: paramSet(
        "LNA2"), background="white", foreground="Blue")
    LNAMenu.add_command(label="Symbolic, Macroscopic", command=lambda: paramSet(
        "LNA3"), background="white", foreground="Blue")
    LNAMenu.add_command(label="Symbolic, f(xo), Macroscopic", command=lambda: paramSet(
        "LNA-xo"), background="white", foreground="Blue")
    LNAMenu.add_command(label="Symbolic, f(ks), Macroscopic", command=lambda: paramSet(
        "LNA-ks"), background="white", foreground="Blue")
    LNAMenu.add_command(label="Symbolic, values, Macroscopic", command=lambda: paramSet(
        "LNA-vs"), background="white", foreground="Blue")
    menubut2.menu.add_cascade(label="Linear Noise Appx.", menu=LNAMenu)

    GilMenu = gui.Menu(frame, tearoff=1)
    GilMenu.add_command(label="Direct method", command=lambda: paramSet(
        "Gillespie_"), background="white", foreground="Blue")
    #GilMenu.add_command ( label="First Rxn Method",command=lambda: print("Not implemented yet"),background="white",foreground="Blue"  )
    #GilMenu.add_command ( label="Next Rxn Method",command=lambda: print("Not implemented yet"),background="white",foreground="Blue"  )
    #GilMenu.add_command ( label="Optimized Direct",command=lambda: print("Not implemented yet"),background="white",foreground="Blue"  )
    menubut2.menu.add_cascade(label="Gillespie", menu=GilMenu)
    menubut2.place(x=95, y=5)

    menubut3 = gui.Menubutton(frame, text="    Analysis    ", activebackground="#f2f20d",
                              activeforeground="red", bg="#00cc00", fg="white" if platform.lower() != "darwin" else "green")
    menubut3.menu = gui.Menu(menubut3, tearoff=1)
    menubut3["menu"] = menubut3.menu
    menubut3.menu.add_command(label="Covariance", command=lambda: analysis_case(
        "cov", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Fano Factor", command=lambda: analysis_case(
        "fanoF", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Cross correlation", command=lambda: analysis_case(
        "corr", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Probability density", command=lambda: analysis_case(
        "pdens1", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Freq. Dist w/r to t", command=lambda: analysis_case(
        "pdens2", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Hist. slice of time", command=lambda: analysis_case(
        "pdens3", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Average of traj.", command=lambda: analysis_case(
        "avetrj", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Phase portrait", command=lambda: analysis_case(
        "phaseP", items), background="white", foreground="Blue")
    menubut3.menu.add_command(label="Plot Data", command=lambda: analysis_case(
        "plotD", items), background="white", foreground="Blue")
    NetLMenu = gui.Menu(frame, tearoff=1)
    NetLMenu.add_command(label="Symbolic, Macroscopic", command=lambda: paramSet(
        "NetLoc1"), background="white", foreground="Blue")
    NetLMenu.add_command(label="Numeric, Macroscopic", command=lambda: paramSet(
        "NetLoc2"), background="white", foreground="Blue")
    menubut3.menu.add_cascade(label="Network Localization", menu=NetLMenu)
    menubut3.place(x=189, y=5)

    frame1 = gui.Frame(frame, height=435, width=972,
                       bg='#8c8c8c', borderwidth=2)
    frame1.place(x=0, y=35, relheight=0.93, relwidth=1.0)
    items = prepare_frame_for_plot(frame1, 972, 435)

    top.mainloop()
