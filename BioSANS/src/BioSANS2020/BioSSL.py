#import sys
import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import pathlib
from BioSANS2020.prepcodes.process import *
from BioSANS2020.cli_functs.ssl_calls import *
import pandas as pd
from BioSANS2020.myglobal import mglobals as globals2
import matplotlib.pyplot as plt
import webbrowser


if 'USERPROFILE' in os.environ:
    cwd = os.path.join(os.environ['USERPROFILE'], "BioSSL_temporary_folder")
elif 'HOME' in os.environ:
    cwd = os.path.join(os.environ['HOME'], "BioSSL_temporary_folder")
else:
    try:
        import tempfile
        cwd = os.path.join(tempfile.gettempdir(), "BioSSL_temporary_folder")
    except BaseException:
        cwd = os.path.join(os.getcwd(), "BioSSL_temporary_folder")

try:
    os.mkdir(cwd, 0o777)
except BaseException:
    pass

globals2.init()
trj = {}


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


def get_input():
    try:
        return raw_input("> ")
    except BaseException:
        return input("> ")


def process_command(command):
    global cwd, trj
    rowc = command.strip().split()
    miter = 1
    tn = 100
    Vol = 1.0
    tsc = 1.5
    tlen = 100
    mult_proc = False
    fout = "temp_traj"
    fileUnit = "molecules"
    plot = True
    mixp = True
    norm = False
    logx = False
    logy = False
    EdataFile = None
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
        fFile = cwd + "/temp.txt"
        f = open(fFile, "w")
        f.write("#REACTIONS\n")
        for rx in rxns:
            f.write(rx + "\n")
        f.write("\n")
        Mname = ""
        f.write("@CONCENTRATION\n")
        conc = ""

        i = i + 1
        while rowc[i] != "using" and i < rlen:
            conc = conc + rowc[i]
            Mname = rowc[i + 2]
            i = i + 1

        for cx in conc.split("&"):
            f.write(cx.replace("=", ",") + "\n")
        f.write("\n")
        f.close()

        i = i + 3
        opts = ""
        while i < rlen:
            opts = opts + rowc[i]
            i = i + 1
        opts = opts.split("&")
        optsv = {}
        for op in opts:
            opr = op.split("=")
            try:
                optsv[opr[0].strip()] = opr[1]
            except BaseException:
                pass

        if 'tn' in optsv:
            tn = float(optsv['tn'])
        if 'tlen' in optsv:
            tlen = int(optsv['tlen'])
        if 'Vol' in optsv:
            Vol = float(optsv['Vol'])
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
            fileUnit = optsv['fileUnit'].lower()
        if 'EdataFile' in optsv:
            EdataFile = optsv['EdataFile']
            EdataFile = os.path.join(cwd, EdataFile)
        if 'topo' in optsv:
            topo = optsv['topo']
            topo = os.path.join(cwd, topo)
            fFile = topo

        Fname = cwd + "/" + fout
        process(
            rfile=fFile,
            miter=miter,
            inMolar=fileUnit,
            Vm=Vol,
            tn=tn,
            delX=tsc,
            normalize=norm,
            logx=logx,
            logy=logy,
            method=Mname,
            tlen=tlen,
            mix_plot=mixp,
            save=True,
            out_fname=Fname,
            plot_show=plot,
            Cinput={},
            vary="",
            mult_proc=mult_proc,
            items=None,
            expDataFile=EdataFile
        )
        # except:
        #print("temp.txt file not found")
    elif rowc[0].lower() == "load":
        sslfile = cwd + "/" + rowc[1]
        f = open(sslfile)
        command = ""
        for row in f:
            if len(row.strip()) > 0:
                if row.strip()[-1] == ";":
                    command = command + row.strip().replace(";", "")
                    process_command(command)
                    command = ""
                else:
                    command = command + row.strip() + " "
        f.close()
    elif rowc[0].lower() == "pwd":
        print(cwd)
    elif rowc[0].lower() == "mkdir":
        try:
            os.mkdir(cwd + "/" + rowc[1], 0o777)
        except BaseException:
            pass
    elif rowc[0].lower() == "ls":
        dirs = []
        try:
            if len(rowc) > 1:
                dirs = os.listdir(cwd + "/" + rowc[1])
            else:
                dirs = os.listdir(cwd)
        except BaseException:
            pass
        for x in dirs:
            if os.path.isdir(x):
                print("directory : " + x)
            else:
                print("file      : " + x)
    elif rowc[0].lower() == "cd":
        try:
            if os.path.isdir(cwd + "/" + rowc[1]):
                abspath = os.path.abspath(cwd + "/" + rowc[1])
                dname = os.path.dirname(abspath)
                os.chdir(dname)
                cwd = str(abspath)
            elif os.path.isdir(rowc[1].strip()):
                abspath = os.path.abspath(rowc[1])
                dname = os.path.dirname(abspath)
                os.chdir(dname)
                cwd = str(abspath)
            else:
                print("No such directory")
        except BaseException:
            print("cannot change dir")
    elif rowc[0].lower() == "read_traj":
        if len(rowc) == 4:
            if rowc[2].lower() == "as":
                name = rowc[3].strip()
                try:
                    sslfile = cwd + "/" + rowc[1]
                    trj[name] = load_data_traj(sslfile)
                except BaseException:
                    try:
                        sslfile = rowc[1]
                        trj[name] = load_data_traj(sslfile)
                    except BaseException:
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
                    print(trj[rowc[1].strip()][eval(rowc[2].strip())])
                except BaseException:
                    print(trj[rowc[1].strip()]
                          [eval("'" + rowc[2].strip() + "'")])
            else:
                print(trj[rowc[1].strip()])
        except BaseException:
            pass
    elif rowc[0].lower() == "plot":
        if len(rowc) == 4:
            trj[rowc[1].strip()].plot(x=rowc[2], y=rowc[3], kind='scatter', s=1)
            plt.show()
        elif len(rowc) > 4:
            qx = rowc[3:]
            plt.xlabel(rowc[2])
            for x in qx:
                plt.scatter(trj[rowc[1]][rowc[2]],
                            trj[rowc[1]][x], label=x, s=1)
            plt.legend()
            plt.show()
        else:
            print("Invalid syntax")
    elif rowc[0].lower() == "calc_covariance":
        try:
            calc_covariance(trj[rowc[1].strip()], int(rowc[2]))
        except BaseException:
            calc_covariance(trj[rowc[1].strip()], 100)
    elif rowc[0].lower() == "prob_density":
        prob_density_calc(trj[rowc[1].strip()], cwd + "/" + rowc[1].strip())
    elif rowc[0].lower() == "prob_density_wtime":
        prob_density_calc_wtime(trj[rowc[1].strip()], cwd +
                                "/" + rowc[1].strip(), "Mname")
    elif rowc[0].lower() == "calc_average":
        try:
            calc_average_conc_at_tend(trj[rowc[1].strip()], int(rowc[2]))
        except BaseException:
            calc_average_conc_at_tend(trj[rowc[1].strip()], 100)
    elif rowc[0].lower() == "length":
        try:
            print(len(trj[rowc[1].strip()]))
        except BaseException:
            pass
    elif rowc[0].lower() == "pdread_traj":
        if len(rowc) == 4:
            if rowc[2].lower() == "as":
                name = rowc[3].strip()
                try:
                    sslfile = cwd + "/" + rowc[1]
                    trj[name] = pd.read_csv(sslfile, delimiter="\t")
                except BaseException:
                    try:
                        sslfile = rowc[1]
                        trj[name] = pd.read_csv(sslfile, delimiter="\t")
                    except BaseException:
                        print("file not found")
            else:
                print("use as to assign content to variable")
                print("read_traj filename as variable")
        else:
            print("invalid syntax")
    elif rowc[0].lower() == "open_pwd":
        show_file_dir(cwd + "/.")


if __name__ == '__main__':
    print("###############################################################")
    print("Welcome to BioSSL commandline interface\n")

    abspath = cwd  # os.path.abspath(cwd)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    cwd = str(abspath).replace("\\", "/")

    row = " "
    command = " "

    while row.strip() != "quit()":
        row = " " + get_input().strip()
        if row[-1] == ";":
            command = command + " " + row.strip().replace(";", "")
            if command.strip() != "":
                process_command(command)
            command = " "
            row = " "
        else:
            command = command + " " + row
