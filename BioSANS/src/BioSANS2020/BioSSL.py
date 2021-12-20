"""

                     This module is the BioSSl module

This module process SSL (structured simulation language) queries.


"""

# import sys
import os
import webbrowser
from sys import platform as PLATFORM
import pandas as pd
import matplotlib.pyplot as plt
# sys.path.append(os.path.abspath("BioSANS2020"))

# import pathlib
from BioSANS2020.prepcodes.process import process
from BioSANS2020.cli_functs.ssl_calls \
    import load_data_traj, calc_average_conc_at_tend, calc_covariance, \
    prob_density_calc_wtime, prob_density_calc  # , \
# calc_covariance_per_traj, calc_covariance_bootsrap, \
# prob_density_calc_tslice
from BioSANS2020.myglobal import mglobals as globals2

if PLATFORM == "win32":
    pass
    # from subprocess import Popen, CREATE_NEW_CONSOLE
elif PLATFORM == "darwin":
    # try:
    # from applescript import tell as my_tell_us
    # except:
    # pass
    from subprocess import call as my_call_us
elif PLATFORM == "linux":
    pass
else:
    pass
    # from subprocess import Popen

if 'USERPROFILE' in os.environ:
    CWD = os.path.join(os.environ['USERPROFILE'], "BioSSL_temporary_folder")
elif 'HOME' in os.environ:
    CWD = os.path.join(os.environ['HOME'], "BioSSL_temporary_folder")
else:
    try:
        import tempfile
        CWD = os.path.join(tempfile.gettempdir(), "BioSSL_temporary_folder")
    except:
        CWD = os.path.join(os.getcwd(), "BioSSL_temporary_folder")

try:
    os.mkdir(CWD, 0o777)
except:
    pass

globals2.init(globals2)
TRAJ = {}


def eval2(to_eval):
    """Evaluate expression"""
    return eval(to_eval)


def show_file_dir(path):
    """This function opens the current working directory

    Args:
        path (str): current working directory
    """
    if PLATFORM == "win32":
        os.startfile(os.path.dirname(path))
    elif PLATFORM == "darwin":
        my_call_us('open', os.path.dirname(path))
    elif PLATFORM == "linux":
        my_call_us('xdg-open', os.path.dirname(path))
    else:
        webbrowser.open(os.path.dirname(path))


def get_input():
    """Thi function get the input commands as user type in console

    Returns:
        str: The input command
    """
    try:
        return input("> ")
    except:
        return raw_input("> ")


def process_command(command):
    """This function handles the processing of commands.

    Args:
        command (str): string of commands
    """
    global CWD
    rowc = command.strip().split()
    miter = 1
    tend = 100
    vol = 1.0
    tsc = 1.5
    tlen = 100
    mult_proc = False
    fout = "temp_traj"
    file_unit = "molecules"
    plot = True
    mixp = True
    norm = False
    logx = False
    logy = False
    edatafile = None
    topo = None
    if rowc[0].lower() == "propagate":
        rlen = len(rowc)
        i = 1
        rxns = ""
        while rowc[i] != "where" and i < rlen:
            rxns = rxns + " " + rowc[i]
            i = i + 1
        rxns = rxns.split("&")
        # try:
        ffile = CWD + "/temp.txt"
        fvar = open(ffile, "w")
        fvar.write("#REACTIONS\n")
        for r_x in rxns:
            fvar.write(r_x + "\n")
        fvar.write("\n")
        mname = ""
        fvar.write("@CONCENTRATION\n")
        conc = ""

        i = i + 1
        while rowc[i] != "using" and i < rlen:
            conc = conc + rowc[i]
            mname = rowc[i + 2]
            i = i + 1

        for c_x in conc.split("&"):
            fvar.write(c_x.replace("=", ",") + "\n")
        fvar.write("\n")
        fvar.close()

        i = i + 3
        opts = ""
        while i < rlen:
            opts = opts + rowc[i]
            i = i + 1
        opts = opts.split("&")
        optsv = {}
        for o_p in opts:
            opr = o_p.split("=")
            try:
                optsv[opr[0].strip()] = opr[1]
            except:
                pass

        if 'tend' in optsv:
            tend = float(optsv['tend'])
        if 'tlen' in optsv:
            tlen = int(optsv['tlen'])
        if 'Vol' in optsv:
            vol = float(optsv['Vol'])
        if 'tsc' in optsv:
            tsc = float(optsv['tsc'])
        if 'miter' in optsv:
            miter = int(optsv['miter'])
        if 'mult_proc' in optsv:
            mult_proc = optsv['mult_proc'].lower() == "true"
        if 'fout' in optsv:
            fout = optsv['fout']
        if 'plot' in optsv:
            plot = optsv['plot'].lower() == "true"
        if 'norm' in optsv:
            norm = optsv['norm'].lower() == "true"
        if 'mixp' in optsv:
            mixp = optsv['mixp'].lower() == "true"
        if 'logx' in optsv:
            norm = optsv['norm'].lower() == "logx"
        if 'logy' in optsv:
            norm = optsv['norm'].lower() == "logy"
        if 'fileUnit' in optsv:
            file_unit = optsv['fileUnit'].lower()
        if 'EdataFile' in optsv:
            edatafile = optsv['EdataFile']
            edatafile = os.path.join(CWD, edatafile)
        if 'topo' in optsv:
            topo = optsv['topo']
            topo = os.path.join(CWD, topo)
            ffile = topo

        fname = CWD + "/" + fout
        process(
            rfile=ffile,
            miter=miter,
            conc_unit=file_unit,
            v_volms=vol,
            tend=tend,
            del_coef=tsc,
            normalize=norm,
            logx=logx,
            logy=logy,
            method=mname,
            tlen=tlen,
            mix_plot=mixp,
            save=True,
            out_fname=fname,
            plot_show=plot,
            c_input={},
            vary="",
            mult_proc=mult_proc,
            items=None,
            exp_data_file=edatafile
        )
        # except:
        #print("temp.txt file not found")
    elif rowc[0].lower() == "load":
        sslfile = CWD + "/" + rowc[1]
        fvar = open(sslfile)
        command = ""
        for row in fvar:
            if len(row.strip()) > 0:
                if row.strip()[-1] == ";":
                    command = command + row.strip().replace(";", "")
                    process_command(command)
                    command = ""
                else:
                    command = command + row.strip() + " "
        fvar.close()
    elif rowc[0].lower() == "pwd":
        print(CWD)
    elif rowc[0].lower() == "mkdir":
        try:
            os.mkdir(CWD + "/" + rowc[1], 0o777)
        except:
            pass
    elif rowc[0].lower() == "ls":
        dirs = []
        try:
            if len(rowc) > 1:
                dirs = os.listdir(CWD + "/" + rowc[1])
            else:
                dirs = os.listdir(CWD)
        except:
            pass
        for xvar in dirs:
            if os.path.isdir(xvar):
                print("directory : " + xvar)
            else:
                print("file      : " + xvar)
    elif rowc[0].lower() == "cd":
        try:
            if os.path.isdir(CWD + "/" + rowc[1]):
                abspath = os.path.abspath(CWD + "/" + rowc[1])
                dname = os.path.dirname(abspath)
                os.chdir(dname)
                CWD = str(abspath)
            elif os.path.isdir(rowc[1].strip()):
                abspath = os.path.abspath(rowc[1])
                dname = os.path.dirname(abspath)
                os.chdir(dname)
                CWD = str(abspath)
            else:
                print("No such directory")
        except:
            print("cannot change dir")
    elif rowc[0].lower() == "read_traj":
        if len(rowc) == 4:
            if rowc[2].lower() == "as":
                name = rowc[3].strip()
                try:
                    sslfile = CWD + "/" + rowc[1]
                    TRAJ[name] = load_data_traj(sslfile)
                except:
                    try:
                        sslfile = rowc[1]
                        TRAJ[name] = load_data_traj(sslfile)
                    except:
                        print("File not found")
            else:
                print("use as to assign content to variable")
                print("read_traj filename as variable")
        else:
            print("Invalid syntax")
    elif rowc[0].lower() == "print":
        try:
            if len(rowc) == 3:
                try:
                    print(TRAJ[rowc[1].strip()][eval2(rowc[2].strip())])
                except:
                    print(TRAJ[rowc[1].strip()]
                          [eval2("'" + rowc[2].strip() + "'")])
            else:
                print(TRAJ[rowc[1].strip()])
        except:
            pass
    elif rowc[0].lower() == "plot":
        if len(rowc) == 4:
            TRAJ[rowc[1].strip()].plot(
                xvar=rowc[2], y=rowc[3], kind='scatter', s=1)
            plt.show()
        elif len(rowc) > 4:
            q_x = rowc[3:]
            plt.xlabel(rowc[2])
            for xvar in q_x:
                plt.scatter(TRAJ[rowc[1]][rowc[2]],
                            TRAJ[rowc[1]][xvar], label=xvar, s=1)
            plt.legend()
            plt.show()
        else:
            print("Invalid syntax")
    elif rowc[0].lower() == "calc_covariance":
        try:
            calc_covariance(TRAJ[rowc[1].strip()], int(rowc[2]))
        except:
            calc_covariance(TRAJ[rowc[1].strip()], 100)
    elif rowc[0].lower() == "prob_density":
        prob_density_calc(TRAJ[rowc[1].strip()], CWD + "/" + rowc[1].strip())
    elif rowc[0].lower() == "prob_density_wtime":
        prob_density_calc_wtime(TRAJ[rowc[1].strip()], CWD +
                                "/" + rowc[1].strip(), "Mname")
    elif rowc[0].lower() == "calc_average":
        try:
            calc_average_conc_at_tend(TRAJ[rowc[1].strip()], int(rowc[2]))
        except:
            calc_average_conc_at_tend(TRAJ[rowc[1].strip()], 100)
    elif rowc[0].lower() == "length":
        try:
            print(len(TRAJ[rowc[1].strip()]))
        except:
            pass
    elif rowc[0].lower() == "pdread_traj":
        if len(rowc) == 4:
            if rowc[2].lower() == "as":
                name = rowc[3].strip()
                try:
                    sslfile = CWD + "/" + rowc[1]
                    TRAJ[name] = pd.read_csv(sslfile, delimiter="\t")
                except:
                    try:
                        sslfile = rowc[1]
                        TRAJ[name] = pd.read_csv(sslfile, delimiter="\t")
                    except:
                        print("file not found")
            else:
                print("use as to assign content to variable")
                print("read_traj filename as variable")
        else:
            print("invalid syntax")
    elif rowc[0].lower() == "open_pwd":
        show_file_dir(CWD + "/.")


if __name__ == '__main__':
    print("###############################################################")
    print("Welcome to BioSSL commandline interface\n")

    ABSPATH = CWD  # os.path.abspath(CWD)
    DNAME = os.path.dirname(ABSPATH)
    os.chdir(DNAME)
    CWD = str(ABSPATH).replace("\\", "/")

    ROW = " "
    COMMAND = " "

    while ROW.strip() != "quit()":
        ROW = " " + get_input().strip()
        if ROW[-1] == ";":
            COMMAND = COMMAND + " " + ROW.strip().replace(";", "")
            if COMMAND.strip() != "":
                process_command(COMMAND)
            COMMAND = " "
            ROW = " "
        else:
            COMMAND = COMMAND + " " + ROW
