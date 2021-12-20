"""

                   This module is the BioSANS module

The purpose of this module  is to provide a graphical user interface and
to facilitate  the transfer of information from BioSANS topology file to
other modules.


The following are the list of functions in this module

1. load_data
2. show_file_dir
3. create_file
4. extract_ode
5. sbml_to_topo2
6. save_file
7. runpy_file
8. run_ssl
9. load_data2
10. tload_data2
11. delete_this
12. canvas_update_widgets
13. load_image
14. eval2
15. dict_trans
16. convert
17. range_trans
18. range_prep
19. mrun_propagation
20. tprocess
21. analysis_case
22. plot_traj_d2
23. param_set


"""

import os
import tkinter as gui
from tkinter import ttk  # , Button, filedialog
import time
import webbrowser
from sys import platform as PLATFORM, executable as sys_executable
from math import ceil as Myceil
from pathlib import Path
from datetime import datetime
import threading
from PIL import Image as Image2, ImageTk
from numpy import array as np_array

# from queue import Queue


from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.myglobal import proc_global
from BioSANS2020.gui_functs.prepare_canvas import prepare_frame_for_plot
from BioSANS2020.prepcodes.process import process
from BioSANS2020.analysis.plotting.plot_traj \
    import plot_traj, plot_traj2
from BioSANS2020.analysis.numeric.transform_data import calc_cross_corr, \
    calc_covariance, fano_factor, prob_density_calc, \
    prob_density_calc2, prob_density_calc3, ave_traj_calc
from BioSANS2020.model.fileconvert.process_sbml \
    import process_sbml as sbml_to_topo
from BioSANS2020.model.ode_parse import ode_extract

import BioSANS2020.model.topology_view as topology_view
from BioSANS2020.model.new_file import new_file

if PLATFORM == "win32":
    from subprocess import Popen, CREATE_NEW_CONSOLE
elif PLATFORM == "darwin":
    from applescript import tell as my_tell_us
    from subprocess import check_output, call as my_call_us
elif PLATFORM == "linux":
    pass
else:
    from subprocess import Popen

try:
    import tempfile
    TEMPORARY_FOLDER = str(tempfile.gettempdir()).replace(
        "\\", "/") + "/BioSANS_temporary_folder"
except BaseException:
    TEMPORARY_FOLDER = "BioSANS_temporary_folder"

globals2.init(globals2)
if __name__ == '__main__':
    proc_global.init(proc_global)

TOP = gui.Tk()
TOP.title("BioSANS 1.0")
TOP.geometry("1005x550")
# TOP.resizable(False, False)

HEADER = gui.Label(TOP, text="BioSANS")
HEADER.configure(
    bg="green",
    fg="white",
    height=1,
    # width = 1005,
    font="Helvetica 18 bold italic"
)
HEADER.pack(fill="x")

FRAME = gui.Frame(TOP)
FRAME.configure(
    bg="light cyan",
    borderwidth=2,
    height=500,
    width=1005
)
FRAME.pack(fill="both", expand=True)

FOOTER = gui.Label(
    TOP, text="Biological Symbolic and Numeric Simulation Algorithms")
FOOTER.configure(
    bg="green",
    fg="white",
    # width = 1005,
    font="Helvetica 10 bold italic",
    anchor='w'
)
FOOTER.pack(fill="x")

FILE_NAME = {"topology": TEMPORARY_FOLDER, "current_folder": None}


def load_data(itups):
    """This fuction reads topology file and display the contents in a
    text area.

    Args:
        itups (tuple): (canvas, scroll_x, scroll_y)
    """
    global FILE_NAME
    file = gui.filedialog.askopenfilename(title="Select file")
    FILE_NAME["topology"] = file
    FILE_NAME["current_folder"] = file
    globals2.TO_CONVERT = file
    if os.path.isfile(file):
        FILE_NAME['last_open'] = topology_view.view_topo(file, itups)


def show_file_dir(path):
    """This function opens the current working directory

    Args:
        path (str): directory path
    """
    global PLATFORM
    if PLATFORM == "win32":
        os.startfile(os.path.dirname(path))
    elif PLATFORM == "darwin":
        my_call_us('open', os.path.dirname(path))
    elif PLATFORM == "linux":
        my_call_us('xdg-open', os.path.dirname(path))
    else:
        webbrowser.open(os.path.dirname(path))


def create_file(itups, ftype):
    """This function creates a temporary file in a temporary directory
    where new topology / files can be placed.

    Args:
        itups (tuple): (canvas, scroll_x, scroll_y)
        ftype (int): 1 for BioSANS topology file, 2 for ODE file
    """
    global FILE_NAME, TEMPORARY_FOLDER
    try:
        os.mkdir(TEMPORARY_FOLDER, 0o777)
    except BaseException:
        for item in Path(TEMPORARY_FOLDER).iterdir():
            if item.is_dir():
                pass
            else:
                item.unlink()

    FILE_NAME['last_open'] = new_file(itups)
    FILE_NAME["topology"] = TEMPORARY_FOLDER + "/temp.txt"
    if ftype == 1:
        FILE_NAME['last_open'].insert('insert', "FUNCTION_DEFINITIONS:\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert(
            'insert', "#REACTIONS, Volume = 1, tend = 100, steps = 100, \
            FileUnit = molar\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "@CONCENTRATION\n")
    elif ftype == 2:
        FILE_NAME['last_open'].insert('insert', "ODE_DECLARATIONS:\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "INI_CONCENTRATIONS:\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "\n")
        FILE_NAME['last_open'].insert('insert', "RATE_CONSTANTS:\n")


def extract_ode(itups):
    """This function extracts topology from ODE file format

    Args:
        itups (tuple): (canvas, scroll_x, scroll_y)
    """
    # global FILE_NAME
    FILE_NAME['last_open'] = ode_extract.odedxdt_to_topo(
        FILE_NAME["topology"], itups)
    FILE_NAME["topology"] = FILE_NAME["topology"] + ".TOP"


def sbml_to_topo2(tocon, itups):
    """This function helps  in the conversion  of SBML files  to BioSANS
    topology files.

    Args:
        tocon (str): sbml file name to convert
        itups (tuple): (canvas, scroll_x, scroll_y)
    """
    FILE_NAME["topology"] = sbml_to_topo(tocon)
    FILE_NAME['last_open'] = topology_view.view_topo(
        FILE_NAME["topology"], itups)


def save_file():
    """This function saves the content of the last opened text area into
    a file.
    """
    # global FILE_NAME
    file = gui.filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    FILE_NAME["topology"] = file.name
    FILE_NAME["current_folder"] = file.name
    if file is None:
        return
    file.write(FILE_NAME['last_open'].get("0.0", "end"))
    file.close()


def runpy_file():
    # global FILE_NAME, PIPE
    with open(FILE_NAME["topology"], "w") as ffvar:
        ffvar.write(FILE_NAME['last_open'].get("0.0", "end"))
        ffvar.write("\ninput('Press enter to exit:')")
    if PLATFORM == "win32":
        Popen([sys_executable, FILE_NAME["topology"]],
              creationflags=CREATE_NEW_CONSOLE)
    elif PLATFORM == "darwin":
        my_tell_us.app("Terminal", 'do script "' +
                       str(sys_executable) + " " + FILE_NAME["topology"] + '"')
    elif PLATFORM == "linux":
        os.system("gnome-terminal -x python3 " + FILE_NAME["topology"])
    else:
        Popen(str(sys_executable) + " " + FILE_NAME["topology"], shell=True)


def run_ssl():
    """This function initiates BioSSL or the BioSABS structured
    simulation language."""
    if PLATFORM == "win32":
        Popen([sys_executable, "-m", "BioSANS2020.BioSSL"],
              creationflags=CREATE_NEW_CONSOLE)
    elif PLATFORM == "darwin":
        pip_show_out = check_output(['pip3', 'show', 'BioSANS2020-efajiculay'])
        pip_show_out = str(pip_show_out).split("\\n")
        install_dir = ""
        for row in pip_show_out:
            line = row.split(":")
            if line[0].strip() == "Location":
                install_dir = (
                    "".join(line[1:])
                    .strip("\\r\\n").replace("c\\", "c:/")
                    .replace("\\", "/").replace("//", "/"))
        if install_dir != "":
            install_dir = str(install_dir)
            my_tell_us.app("Terminal", 'do script "' + str(sys_executable) +
                           " " + install_dir + "/BioSANS2020/BioSSL.py" + '"')
    elif PLATFORM == "linux":
        os.system("gnome-terminal -x python3 " +
                  os.path.join(os.getcwd(), "BioSSL.py"))
    else:
        Popen(str(sys_executable) + " " +
              os.path.join(os.getcwd(), "BioSSL.py"), shell=True)


def load_data2(plot=False):
    """This function loads data for numerical processing or for plotting

    Args:
        plot (bool, optional): if True, loaded trajectory is plotted.
        Defaults to False.
    """
    # t_o = time.time()
    global CURRENT_DATA
    file = gui.filedialog.askopenfilename(title="Select file")
    FILE_NAME["trajectory"] = file
    FILE_NAME["current_folder"] = file
    with open(FILE_NAME["trajectory"], "r") as fvar:
        data = []
        ddvar = []
        row1 = str(fvar.readline()).strip()
        slabels = row1.split("\t")[1:]
        for row in fvar:
            cols = [float(xvar) for xvar in row.split("\t")]
            if cols[0] == 0.0 and ddvar:
                data.append(np_array(ddvar))
                ddvar = []
            ddvar.append(cols)
        data.append(np_array(ddvar))
    if plot:
        plot_traj(data, slabels, ITUPS, globals2.PLOTTED, mix_plot=True,
                  logx=False, logy=False, normalize=False)
    # print(data[0], "\n\n")
    # print(data[1], "\n\n")
    # print(data[2], "\n")
    CURRENT_DATA = (data, slabels)
    gui.messagebox.showinfo("showinfo", "Trajectory loaded succesfully")
    # print(time.time()-t_o)


def tload_data2(plot=False):
    """This function starts a thread that handles loading data.

    Args:
        plot (bool, optional): if True, loaded trajectory is plotted.
        Defaults to False.
    """
    if __name__ == '__main__':
        tvar = threading.Thread(target=load_data2, args=(plot,), daemon=False)
        tvar.start()


def delete_this(frame, canvas):
    canvas.delete(frame)
    canvas_update_widgets(None, canvas)


def canvas_update_widgets(_, canvas):
    # R = canvas._root().winfo_height()
    # root = canvas._root()
    # width = root.winfo_screenwidth()
    # height = root.winfo_screenheight()
    # print(canvas.winfo_children())
    cheigh = canvas.winfo_height()
    cwidth = canvas.winfo_width()
    obj_can = canvas.find_all()
    objLen = len(obj_can)
    ind = 0
    for xvar in obj_can:
        canvas.itemconfig(xvar, height=cheigh, width=cwidth - 5)
        # x_1, y_1, x_2, y_2 = canvas.bbox( xvar)
        _, y_1, _, _ = canvas.bbox(xvar)
        canvas.move(xvar, 0, (objLen - ind) * cheigh - y_1 + 3)
        ind = ind + 1
    canvas.configure(scrollregion=canvas.bbox("all"))
    return "break"


def load_image(wdata=False):
    t_o = time.time()
    global CURRENT_DATA
    canvas, scroll_x, scroll_y = ITUPS
    file = gui.filedialog.askopenfilename(title="Select file")
    FILE_NAME["current_folder"] = file
    load = Image2.open(file)
    render = ImageTk.PhotoImage(load)
    img = gui.Label(canvas, image=render)
    img.image = render

    fframe = canvas.create_window(
        0, 426 * globals2.PLOT_I, anchor='nw', window=img)
    but_b = gui.Button(img, text=" X ", fg='red', highlightcolor='blue',
                       bg='white', height=1, relief='raised',
                       command=lambda: delete_this(fframe, canvas))
    but_b.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")

    canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
    canvas.configure(scrollregion=canvas.bbox("all"))
    canvas.bind("<Configure>", lambda e: canvas_update_widgets(e, canvas))
    globals2.PLOT_I = globals2.PLOT_I + 1

    if wdata:
        file = str(file).replace("jpg", "dat").replace("png", "dat")
        with open(file, "r") as ffile:
            data = []
            ddvar = []
            row1 = str(ffile.readline()).strip()
            slabels = row1.split("\t")[1:]
            for row in ffile:
                cols = [float(xvar) for xvar in row.split("\t")]
                if cols[0] == 0.0 and ddvar:
                    data.append(np_array(ddvar))
                    ddvar = []
                ddvar.append(cols)
            data.append(np_array(ddvar))
        CURRENT_DATA = (data, slabels)
    print(time.time() - t_o)
    canvas_update_widgets(None, canvas)


def eval2(xvar):
    try:
        return eval(xvar)
    except BaseException:
        try:
            par = str(xvar).lower().capitalize()
            return eval(par)
        except BaseException:
            return str(xvar)


CURRENT_DATA = None


def dict_trans(x_1):
    x_2 = x_1.split(",")
    x_3 = {}
    try:
        for xvar in x_2:
            rvar = xvar.split("=")
            x_3[rvar[0].strip()] = float(rvar[1])
    except BaseException:
        pass
    return x_3


def convert(xvar, con):
    try:
        return con(xvar)
    except BaseException:
        return xvar


def range_trans(x_1):
    x_3 = []
    x_2 = x_1.split(",")
    if len(x_2) > 1:
        x_3 = x_2[1:]
        x_3.append(x_2[0])
    return x_3


def range_prep(x_1):
    x_3 = []
    x_2 = x_1.split(",")
    if len(x_2) > 1:
        x_2 = [convert(x_2[xvar], float) for xvar in range(len(x_2))]
        c_c = x_2[0].lower().split("f")
        r_2 = 0
        if len(c_c) < 2:
            c_c = x_2[0].lower().split("b")
            r_2 = 1
        r_1 = int(c_c[1]) - 1
        c_c = ",".join([str(xvar) for xvar in x_2[1:]])
        x_3 = list(eval(c_c))
        return [[r_1, r_2], x_3]
    return x_1


def mrun_propagation(par, entry_list, defs2):
    global CURRENT_DATA
    try:
        del FILE_NAME["trajectory"]
    except BaseException:
        pass

    defs = []
    for i in range(len(entry_list)):
        try:
            val = eval2(entry_list[i].get())
            defs.append(val)
        except BaseException:
            val = eval2(defs2[i].get())
            defs.append(val)
    # defs[15] = dict_trans(entry_list[15].get())
    defs[16] = range_trans(entry_list[16].get())
    defs[17] = range_prep(entry_list[17].get())
    if not defs[9] in [
            "Analyt", "SAnalyt", "Analyt-ftx", "SAnalyt-ftk", "k_est1",
            "k_est2", "k_est3", "k_est4", "k_est5", "k_est6", "k_est7",
            "k_est8", "k_est9", "k_est10", "k_est11", "NetLoc1", "NetLoc2"]:
        with open(defs[13] + "_" + defs[9] + "_params.dat", "w") as fvar:
            fvar.write("\n".join([str(xvar) for xvar in defs]))
    if defs[9] in [
            "k_est1", "k_est2", "k_est3", "k_est4", "k_est6", "k_est7",
            "k_est8", "k_est9", "k_est10", "k_est11", "LNA2", "LNA3", "Analyt",
            "Analyt-ftx", "SAnalyt", "SAnalyt-ftk", "Analyt2", "topoTosbml",
            "topoTosbml2", "topoTosbml3", "LNA-vs", "LNA-ks", "LNA-xo",
            "NetLoc1", "NetLoc2"]:
        par.destroy()
    defs.append(ITUPS)
    tprocess(defs)


SUPER_THREAD_RUN = None


def tprocess(defs):
    global SUPER_THREAD_RUN
    if defs[9] == "k_est5":
        process(*defs)
    else:
        if __name__ == '__main__':
            SUPER_THREAD_RUN = 1
            tvar = threading.Thread(
                target=lambda: process(
                    *defs), daemon=False)
            tvar.start()


def analysis_case(ana_case, itups):
    global FILE_NAME
    if "trajectory" not in FILE_NAME:
        gui.messagebox.showinfo(
            "Trajectory not loaded yet", "Please load the trajectory. BioSANS \
            save it into a file during your last run.")
        try:
            load_data2(False)
        except BaseException:
            try:
                del FILE_NAME["trajectory"]
            except BaseException:
                pass
            return None
        # data, slabels = CURRENT_DATA

    if ana_case == "cov":
        return calc_covariance(CURRENT_DATA, itups)
    if ana_case == "fanoF":
        return fano_factor(CURRENT_DATA, itups)
    if ana_case == "corr":
        return calc_cross_corr(CURRENT_DATA, itups)
    if ana_case == "pdens1":
        return prob_density_calc(CURRENT_DATA, itups)
    if ana_case == "pdens2":
        return prob_density_calc2(CURRENT_DATA, itups)
    if ana_case == "pdens3":
        return prob_density_calc3(CURRENT_DATA, itups)
    if ana_case == "avetrj":
        return ave_traj_calc(CURRENT_DATA, itups)
    if ana_case == "phaseP":
        return plot_traj_d(CURRENT_DATA, itups)
    if ana_case == "plotD":
        return plot_traj_d2(CURRENT_DATA, itups)
    return None


def plot_traj_d(current_data, itups):
    # global SUPER_THREAD_RUN
    if SUPER_THREAD_RUN == 1:
        gui.messagebox.showinfo(
            "Warning: Unsafe thread",
            "If you want to have a 3D Phase portrait, restart BioSANS and \
            load trajectory. For 2D phase portrait, just continue")
    try:
        data, slabels = current_data
        par = gui.Toplevel()
        par.resizable(False, False)
        par.wm_title("Plot settings")

        entry_lst = [
            gui.Label(par, text="choose x-axis", fg="blue"),
            gui.Label(par, text="choose y-axis", fg="blue"),
            gui.Label(par, text="choose z-axis", fg="blue"),
            gui.Label(par, text="choose step-range", fg="blue")
        ]
        for i in range(4):
            entry_lst[i].grid(row=i, column=0, sticky=gui.W, pady=2)
        optsel = [[par, gui.StringVar(), "None", "time"] +
                  slabels for i in range(3)]
        opt_sel_var = [optsel[i][1] for i in range(3)]
        for i in range(3):
            opt_sel_var[i].set("None")
        entry_list = [gui.OptionMenu(*optsel[i]) for i in range(3)]
        for i in range(3):
            entry_list[i].config(width=14)
        entry_list.append(gui.Entry(par, bd=5))
        for i in range(4):
            entry_list[i].grid(row=i, column=1, sticky=gui.W, pady=2)
        entry_list[-1].insert(gui.END, "0:-1")

        but_b1 = ttk.Button(par, text="PLOT", command=lambda: plot_traj2(
            data, slabels, itups, globals2.PLOTTED, logx=False, logy=False,
            normalize=False, xlabel=opt_sel_var[0].get(),
            ylabel=opt_sel_var[1].get(), zlabel=opt_sel_var[2].get(),
            trange=entry_list[-1].get()))
        but_b1.grid(row=5, column=0, sticky=gui.W, pady=2)
    except BaseException:
        gui.messagebox.showinfo(
            "Trajectory not loaded yet",
            "Please load the trajectory. BioSANS save it into a file during \
            your last run.")


def get_checked(el_1, slabels):
    check_si = []
    for i, _ in enumerate(el_1):
        key = el_1[i].get()
        if key != "0":
            check_si.append(slabels.index(key))
    return check_si


def plot_traj_d2(current_data, itups):
    try:
        data, slabels = current_data
        pard = gui.Toplevel()
        pard.resizable(True, True)
        pard.wm_title("Plot species")
        pard.maxsize(width=300, height=700)
        canvas, scroll_x, scroll_y = prepare_frame_for_plot(
            pard, width=200, height=200)
        par = gui.Frame(canvas, width=200, height=int(len(slabels) * 150))
        par.pack(side="left", fill="both", expand=True)
        pard.configure(
            bg="light blue",
            borderwidth=2,
        )
        canvas.configure(
            bg="light blue",
            borderwidth=2,
        )

        el_s = len(slabels)
        el_1 = [gui.StringVar() for i in range(el_s)]
        for i in range(el_s):
            el_1[i].set("0")
        entry_lst = [
            gui.Checkbutton(
                par, text=slabels[i], variable=el_1[i], onvalue=slabels[i],
                offvalue="0") for i in range(el_s)]

        for i in range(Myceil(el_s / 2)):
            entry_lst[i].grid(row=i, column=0, sticky=gui.W, pady=2)
        for i in range(Myceil(el_s / 2), el_s):
            entry_lst[i].grid(
                row=i - Myceil(el_s / 2), column=1, sticky=gui.W, pady=2)

        canvas.create_window(0, 0, anchor='n', window=par)
        canvas.configure(yscrollcommand=scroll_y.set,
                         xscrollcommand=scroll_x.set)
        canvas.configure(scrollregion=canvas.bbox("all"))

        but_b1 = ttk.Button(
            pard, text="PLOT",
            command=lambda: plot_traj(
                data, slabels, itups, globals2.PLOTTED, mix_plot=True,
                logx=False, logy=False, normalize=False,
                si_ticked=get_checked(el_1, slabels)))
        but_b1.pack(side="bottom", fill="x")
    except BaseException:
        gui.messagebox.showinfo(
            "Trajectory not loaded yet",
            "Please load the trajectory. BioSANS save it into a file \
                during your last run.")


def param_set(method):
    global FILE_NAME, ITUPS
    with open(FILE_NAME["topology"], "w") as ffvar:
        ffvar.write(FILE_NAME['last_open'].get("0.0", "end"))

    path = Path(FILE_NAME["topology"])
    s_s = str(FILE_NAME["topology"]).split("/")
    s_s = s_s[-1] if len(s_s) > 1 else ""
    name = str(path.parent) + "/" + s_s + "_" + \
        datetime.now().strftime("%Y%m%d_%H%M%S")
    par = gui.Toplevel()
    par.resizable(False, False)
    par.wm_title("Parameter setting")
    opts = [
        "File name",
        "Number of iteration :",
        "File Units? :",
        "Volume (L) :",
        "end time (tend) :",
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
    defs = [
        FILE_NAME["topology"], 1, gui.StringVar(), 1.0, 100, 1.5, False, False,
        False, method, 1000, True, True, name, True, "time (sec)", "", "",
        False, False]
    defs[2].set('molecules')

    topfile = open(FILE_NAME["topology"], "r")
    for row in topfile:
        if row[0] == "#":
            g_g = row.split(",")[1:]
            for xvar in g_g:
                x_x = [g.strip() for g in xvar.split("=")]
                if x_x[0] == "Volume":
                    defs[3] = x_x[1]
                elif x_x[0] == "tend":
                    defs[4] = x_x[1]
                elif x_x[0] == "FileUnit":
                    defs[2].set(x_x[1])
                elif x_x[0] == "logx":
                    defs[7] = x_x[1]
                elif x_x[0] == "Normalized":
                    defs[6] = x_x[1]
                elif x_x[0] == "steps":
                    defs[10] = x_x[1]

    oplen = len(opts)
    entry_lst = [gui.Label(par, text=opts[i], fg="blue") for i in range(oplen)]
    for i in range(10):
        entry_lst[i].grid(row=i, column=0, sticky=gui.W, pady=2)
    for i in range(10):
        entry_lst[10 + i].grid(row=i, column=2, sticky=gui.W, pady=2)\

    entry_list = [gui.Entry(par, bd=5) if i != 2 else gui.OptionMenu(
        par, defs[2], 'molecules', 'molar', 'moles') for i in range(oplen)]
    for i in range(10):
        entry_list[i].grid(row=i, column=1, sticky=gui.W, pady=2)
    for i in range(10):
        entry_list[10 + i].grid(row=i, column=3, sticky=gui.W, pady=2)
    for i in range(oplen):
        if i != 2:
            entry_list[i].insert(gui.END, str(defs[i]))
    entry_list[9].configure(state="disable")
    entry_list[2].config(width=14)
    if method == "ODE":
        entry_list[5].configure(state="disable")
    if method == "ODE2":
        entry_list[5].configure(state="disable")
    elif method == "Gillespie_":
        entry_list[5].configure(state="disable")
    elif method == "CLE":
        entry_list[5].delete(0, gui.END)
        entry_list[5].insert(gui.END, str(10))

    but_b1 = ttk.Button(par, text="RUN",
                        command=lambda: mrun_propagation(par, entry_list,
                                                         defs))
    but_b1.grid(row=oplen, column=0, sticky=gui.W, pady=2)
    if method in ["k_est1", "k_est2", "k_est3", "k_est4", "k_est6", "k_est7",
                  "k_est8", "k_est9", "k_est10", "k_est11", "LNA2", "LNA3",
                  "LNA-vs", "LNA-ks", "LNA-xo", "NetLoc1", "NetLoc2", "Analyt",
                  "SAnalyt-ftk", "SAnalyt", "Analyt-ftx", "Analyt2",
                  "topoTosbml", "topoTosbml2", "topoTosbml3"]:
        but_b1.invoke()


if __name__ == "__main__":

    MENUBUT1 = gui.Menubutton(
        FRAME, text=" File/Model ", activebackground="#f2f20d",
        activeforeground="red", bg="#00cc00",
        fg="white" if PLATFORM.lower() != "darwin" else "green")
    MENUBUT1.menu = gui.Menu(MENUBUT1, tearoff=1)
    MENUBUT1["menu"] = MENUBUT1.menu
    LOADMENU = gui.Menu(FRAME, tearoff=1)
    LOADMENU.add_command(
        label="Topology/File", command=lambda: load_data(ITUPS),
        background="white", foreground="Blue")
    LOADMENU.add_command(
        label="Trajectory file",
        command=tload_data2, background="white", foreground="Blue")
    LOADMENU.add_command(
        label="Traj. w/ plot", command=lambda: tload_data2(True),
        background="white", foreground="Blue")
    LOADMENU.add_command(
        label="Image of plot", command=load_image,
        background="white", foreground="Blue")
    LOADMENU.add_command(
        label="Image w/ data", command=lambda: load_image(True),
        background="white", foreground="Blue")
    LOADMENU.add_command(
        label="Current folder",
        command=lambda: show_file_dir(FILE_NAME["current_folder"]),
        background="white", foreground="Blue")
    MENUBUT1.menu.add_cascade(label="Open", menu=LOADMENU)
    NEWFMENU = gui.Menu(FRAME, tearoff=1)
    NEWFMENU.add_command(label="Blank File", command=lambda: create_file(
        ITUPS, 0), background="white", foreground="Blue")
    NEWFMENU.add_command(label="Topo File", command=lambda: create_file(
        ITUPS, 1), background="white", foreground="Blue")
    NEWFMENU.add_command(label="ODE File", command=lambda: create_file(
        ITUPS, 2), background="white", foreground="Blue")
    MENUBUT1.menu.add_cascade(label="New File", menu=NEWFMENU)
    MENUBUT1.menu.add_command(label="Save File", command=lambda: save_file())
    MENUBUT1.menu.add_command(
        label="Run File.py", command=lambda: runpy_file())
    MENUBUT1.menu.add_command(label="Run SSL", command=lambda: run_ssl())
    CONVMENU = gui.Menu(FRAME, tearoff=1)
    CONVMENU.add_command(label="SBML to Topo", command=lambda: sbml_to_topo2(
        globals2.TO_CONVERT, ITUPS), background="white", foreground="Blue")
    CONVMENU.add_command(label="ODE to Topo", command=lambda: extract_ode(
        ITUPS), background="white", foreground="Blue")
    TOPSBML = gui.Menu(FRAME, tearoff=1)
    TOPSBML.add_command(label="molecules", command=lambda: param_set(
        "topoTosbml"), background="white", foreground="Blue")
    TOPSBML.add_command(label="molar", command=lambda: param_set(
        "topoTosbml2"), background="white", foreground="Blue")
    TOPSBML.add_command(label="no unit", command=lambda: param_set(
        "topoTosbml3"), background="white", foreground="Blue")
    CONVMENU.add_cascade(label="Topo to SBML", menu=TOPSBML)
    MENUBUT1.menu.add_cascade(label="Convert model", menu=CONVMENU)
    PARESMENU = gui.Menu(FRAME, tearoff=1)
    PARESMENU.add_command(
        label="Nelder-Mead (NM), Macroscopic",
        command=lambda: param_set("k_est6"),
        background="white", foreground="Blue")
    PARESMENU.add_command(
        label="Nelder-Mead (NM), Microscopic",
        command=lambda: param_set("k_est7"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="Powell, Macroscopic", command=lambda: param_set("k_est8"),
        background="white", foreground="Blue")
    PARESMENU.add_command(
        label="Powell, Microscopic", command=lambda: param_set("k_est9"),
        background="white", foreground="Blue")
    PARESMENU.add_command(
        label="L-BFGS-B, Macroscopic",
        command=lambda: param_set("k_est10"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="L-BFGS-B, Microscopic",
        command=lambda: param_set("k_est11"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="NM-Diff. Evol., Macroscopic",
        command=lambda: param_set("k_est3"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="NM-Diff. Evol., Microscopic",
        command=lambda: param_set("k_est4"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="Parameter slider/scanner",
        command=lambda: param_set("k_est5"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="MCEM, Macroscopic",
        command=lambda: param_set("k_est1"), background="white",
        foreground="Blue")
    PARESMENU.add_command(
        label="MCEM, Microscopic",
        command=lambda: param_set("k_est2"), background="white",
        foreground="Blue")
    MENUBUT1.menu.add_cascade(label="Estimate Params", menu=PARESMENU)
    MENUBUT1.place(x=2, y=5)

    MENUBUT2 = gui.Menubutton(
        FRAME, text="Propagation", activebackground="#f2f20d",
        activeforeground="red", bg="#00cc00",
        fg="white" if PLATFORM.lower() != "darwin" else "green")
    MENUBUT2.menu = gui.Menu(MENUBUT2, tearoff=1)
    MENUBUT2["menu"] = MENUBUT2.menu
    ANALMENU = gui.Menu(FRAME, tearoff=1)
    ANALMENU.add_command(
        label="Pure Symbolic :f(t,xo,k)",
        command=lambda: param_set("Analyt"), background="white",
        foreground="Blue")
    ANALMENU.add_command(
        label="Semi-Symbolic :f(t)", command=lambda: param_set("SAnalyt"),
        background="white",
        foreground="Blue")
    ANALMENU.add_command(
        label="Semi-Symbolic :f(t,xo)",
        command=lambda: param_set("Analyt-ftx"), background="white",
        foreground="Blue")
    ANALMENU.add_command(label="Semi-Symbolic :f(t,k)", command=lambda: param_set(
        "SAnalyt-ftk"), background="white", foreground="Blue")
    ANALMENU.add_command(label="For wxmaxima", command=lambda: param_set(
        "Analyt2"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="Analytical soln.", menu=ANALMENU)

    ODEMENU = gui.Menu(FRAME, tearoff=1)
    ODEMENU.add_command(label="Molecules(micro)", command=lambda: param_set(
        "ODE-1"), background="white", foreground="Blue")
    ODEMENU.add_command(label="Molar(macro)", command=lambda: param_set(
        "ODE-2"), background="white", foreground="Blue")
    ODEMENU.add_command(label="Mole(macro)", command=lambda: param_set(
        "ODE-3"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="ODE int", menu=ODEMENU)

    RUNGEK4 = gui.Menu(FRAME, tearoff=1)
    RUNGEK4.add_command(label="Molecules(micro)", command=lambda: param_set(
        "rk4-1"), background="white", foreground="Blue")
    RUNGEK4.add_command(label="Molar(macro)", command=lambda: param_set(
        "rk4-2"), background="white", foreground="Blue")
    RUNGEK4.add_command(label="Mole(macro)", command=lambda: param_set(
        "rk4-3"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="RK4-fix-interval", menu=RUNGEK4)

    RUNGEK4A = gui.Menu(FRAME, tearoff=1)
    RUNGEK4A.add_command(label="Molecules(micro)", command=lambda: param_set(
        "rk4-1a"), background="white", foreground="Blue")
    RUNGEK4A.add_command(label="Molar(macro)", command=lambda: param_set(
        "rk4-2a"), background="white", foreground="Blue")
    RUNGEK4A.add_command(label="Mole(macro)", command=lambda: param_set(
        "rk4-3a"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="RK4-tau-adaptive", menu=RUNGEK4A)

    EULRTAU = gui.Menu(FRAME, tearoff=1)
    EULRTAU.add_command(label="Molecules(micro)", command=lambda: param_set(
        "Euler-1"), background="white", foreground="Blue")
    EULRTAU.add_command(label="Molar(macro)", command=lambda: param_set(
        "Euler-2"), background="white", foreground="Blue")
    EULRTAU.add_command(label="Mole(macro)", command=lambda: param_set(
        "Euler-3"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="Euler (tau-adaptive-1)", menu=EULRTAU)

    EULRTAU2 = gui.Menu(FRAME, tearoff=1)
    EULRTAU2.add_command(label="Molecules(micro)", command=lambda: param_set(
        "Euler2-1"), background="white", foreground="Blue")
    EULRTAU2.add_command(label="Molar(macro)", command=lambda: param_set(
        "Euler2-2"), background="white", foreground="Blue")
    EULRTAU2.add_command(label="Mole(macro)", command=lambda: param_set(
        "Euler2-3"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="Euler (tau-adaptive-2)", menu=EULRTAU2)

    CLETAUA = gui.Menu(FRAME, tearoff=1)
    CLETAUA.add_command(label="Molecules(micro)", command=lambda: param_set(
        "CLE"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="CLE (tau-adaptive)", menu=CLETAUA)

    CLETAUA2 = gui.Menu(FRAME, tearoff=1)
    CLETAUA2.add_command(label="Molecules(micro)", command=lambda: param_set(
        "CLE2"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="CLE (cle-fixIntvl)", menu=CLETAUA2)

    TAULMENU = gui.Menu(FRAME, tearoff=1)
    TAULMENU.add_command(label="Tau-leapingV1-micro", command=lambda: param_set(
        "Tau-leaping"), background="white", foreground="Blue")
    TAULMENU.add_command(label="Tau-leapingV2-micro", command=lambda: param_set(
        "Tau-leaping2"), background="white", foreground="Blue")
    TAULMENU.add_command(label="Sim-TauLeap-micro", command=lambda: param_set(
        "Sim-TauLeap"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="Tau-leaping-micro", menu=TAULMENU)

    LNAMENU = gui.Menu(FRAME, tearoff=1)
    LNAMENU.add_command(label="COV-time-dependent, Macroscopic",
                        command=lambda: param_set("LNA(t)"), background="white", foreground="Blue")
    LNAMENU.add_command(label="FF-time-dependent, Macroscopic",
                        command=lambda: param_set("LNA2(t)"), background="white", foreground="Blue")
    LNAMENU.add_command(label="Numeric, values", command=lambda: param_set(
        "LNA"), background="white", foreground="Blue")
    LNAMENU.add_command(label="Symbolic, Microscopic", command=lambda: param_set(
        "LNA2"), background="white", foreground="Blue")
    LNAMENU.add_command(label="Symbolic, Macroscopic", command=lambda: param_set(
        "LNA3"), background="white", foreground="Blue")
    LNAMENU.add_command(label="Symbolic, f(xo), Macroscopic", command=lambda: param_set(
        "LNA-xo"), background="white", foreground="Blue")
    LNAMENU.add_command(label="Symbolic, f(ks), Macroscopic", command=lambda: param_set(
        "LNA-ks"), background="white", foreground="Blue")
    LNAMENU.add_command(label="Symbolic, values, Macroscopic", command=lambda: param_set(
        "LNA-vs"), background="white", foreground="Blue")
    MENUBUT2.menu.add_cascade(label="Linear Noise Appx.", menu=LNAMENU)

    GILMENU = gui.Menu(FRAME, tearoff=1)
    GILMENU.add_command(label="Direct method", command=lambda: param_set(
        "Gillespie_"), background="white", foreground="Blue")
    # GILMENU.add_command ( label="First Rxn Method",command=lambda:
    # print("Not implemented yet"),background="white",foreground="Blue"  )
    # GILMENU.add_command ( label="Next Rxn Method",command=lambda:
    # print("Not implemented yet"),background="white",foreground="Blue"  )
    # GILMENU.add_command ( label="Optimized Direct",command=lambda:
    # print("Not implemented yet"),background="white",foreground="Blue"  )
    MENUBUT2.menu.add_cascade(label="Gillespie", menu=GILMENU)
    MENUBUT2.place(x=95, y=5)

    MENUBUT3 = gui.Menubutton(
        FRAME, text="    Analysis    ", activebackground="#f2f20d",
        activeforeground="red", bg="#00cc00",
        fg="white" if PLATFORM.lower() != "darwin" else "green")
    MENUBUT3.menu = gui.Menu(MENUBUT3, tearoff=1)
    MENUBUT3["menu"] = MENUBUT3.menu
    MENUBUT3.menu.add_command(
        label="Covariance", command=lambda: analysis_case("cov", ITUPS),
        background="white", foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Fano Factor", command=lambda: analysis_case("fanoF", ITUPS),
        background="white", foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Cross correlation",
        command=lambda: analysis_case("corr", ITUPS), background="white",
        foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Probability density",
        command=lambda: analysis_case("pdens1", ITUPS), background="white",
        foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Freq. Dist w/r to t",
        command=lambda: analysis_case("pdens2", ITUPS), background="white",
        foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Hist. slice of time",
        command=lambda: analysis_case("pdens3", ITUPS), background="white",
        foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Average of traj.",
        command=lambda: analysis_case("avetrj", ITUPS), background="white",
        foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Phase portrait",
        command=lambda: analysis_case("phaseP", ITUPS), background="white",
        foreground="Blue")
    MENUBUT3.menu.add_command(
        label="Plot Data", command=lambda: analysis_case("plotD", ITUPS),
        background="white", foreground="Blue")
    NETLMENU = gui.Menu(FRAME, tearoff=1)
    NETLMENU.add_command(
        label="Symbolic, Macroscopic",
        command=lambda: param_set("NetLoc1"), background="white",
        foreground="Blue")
    NETLMENU.add_command(
        label="Numeric, Macroscopic", command=lambda: param_set("NetLoc2"),
        background="white", foreground="Blue")
    MENUBUT3.menu.add_cascade(label="Network Localization", menu=NETLMENU)
    MENUBUT3.place(x=189, y=5)

    FRAME1 = gui.Frame(FRAME, height=435, width=972,
                       bg='#8c8c8c', borderwidth=2)
    FRAME1.place(x=0, y=35, relheight=0.93, relwidth=1.0)
    ITUPS = prepare_frame_for_plot(FRAME1, 972, 435)

    TOP.mainloop()
