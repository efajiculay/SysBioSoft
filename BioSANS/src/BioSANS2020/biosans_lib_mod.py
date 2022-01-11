"""

          This module handles serving BioSANS as a library

The following are the list of methods that can be used in this module

Stochastic (refer to section 10.2.4)

1.    "CLE"            - Molecules(micro), tau-adaptive
2.    "CLE2"           - Molecules(micro), cle-fixIntvl
3.    "Gillespie_"     - Molecules(micro), Direct method
4.    "Tau-leaping"    - Molecules(micro),
                        Not swapping with Gillespie
5.    "Tau-leaping2"   - Molecules(micro),
                        Swapping with Gillespie
6.    "Sim-TauLeap"    - Molecules(micro), Simplified,
                        Swapping with Gillespie

Deterministic (refer to section 10.2.1)

7.    "Euler-1"        - Molecules(micro), tau-adaptive-1
8.    "Euler-2"        - Molar (macro), tau-adaptive-1
9.    "Euler-3"        - Mole (macro), tau-adaptive-1
10.    "Euler2-1"         - Molecules(micro), tau-adaptive-2
11.    "Euler2-2"       - Molar (macro), tau-adaptive-2
12.    "Euler2-3"       - Mole (macro), tau-adaptive-2
13.    "ODE-1"          - Molecules(micro),
                        using ode_int from scipy
14.    "ODE-2"          - Molar(macro),
                        using ode_int from scipy
15.    "ODE-3"          - Mole(macro), using ode_int from scipy
16.    "rk4-1"          - Molecules(micro), fix-interval
17.    "rk4-2"          - Molar(macro), fix-interval
18.    "rk4-3"          - Mole(macro), fix-interval
19.    "rk4-1a"         - Molecules(micro), tau-adaptive
20.    "rk4-2a"         - Molar(macro), tau-adaptive
21.    "rk4-3a"         - Mole(macro), tau-adaptive

Linear Noise Approximation (refer to 10.1.2 & 10.2.2)

22.    "LNA"             - Numeric, values
23.    "LNA-vs"          - Symbolic, values, Macroscopic
24.    "LNA-ks"          - Symbolic, f(ks), Macroscopic
25.    "LNA-xo"          - Symbolic, f(xo), Macroscopic
26.    "LNA2"            - Symbolic, f(xo,ks), Microscopic
27.    "LNA3"            - Symbolic, f(xo,ks), Macroscopic
28.    "LNA(t)"          - COV-time-dependent, Macroscopic
29.    "LNA2(t)"         - FF-time-dependent, Macroscopic

Network Localization (refer to 10.1.3)

30.    "NetLoc1"         - Symbolic, Macroscopic
31.    "NetLoc2"         - Numeric, Macroscopic

Parameter estimation (refer to 10.2.3)

32.    "k_est1"          - MCEM, Macroscopic
33.    "k_est2"          - MCEM, Microscopic
34.    "k_est3"          - NM-Diff. Evol., Macroscopic
35.    "k_est4"          - NM-Diff. Evol., Microscopic
36.    "k_est5"          - Parameter slider/scanner
37.    "k_est6"          - Nelder-Mead (NM), Macroscopic
38.    "k_est7"          - Nelder-Mead (NM), Microscopic
39.    "k_est8"          - Powell, Macroscopic
40.    "k_est9"          - Powell, Microscopic
41.    "k_est10"         - L-BFGS-B, Macroscopic
42.    "k_est11"         - L-BFGS-B, Microscopic

Symbolic/Analytical expression of species (refer to 10.1.1)

43.    "Analyt"          - Pure Symbolic :f(t,xo,k)
44.    "Analyt-ftx"      - Semi-Symbolic :f(t,xo)
45.    "SAnalyt"         - Semi-Symbolic :f(t)
46.    "SAnalyt-ftk"     - Semi-Symbolic :f(t,k)
47.    "Analyt2"         - Creates commands for wxmaxima

"""


import time
import numpy as np
from os import remove as delete_file
from BioSANS2020.prepcodes.process import process
from BioSANS2020.model.fileconvert.process_sbml \
    import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
from BioSANS2020.math_functs.sbml_math import SBML_FUNCT_DICT
from BioSANS2020.propagation.propensity \
    import propensity_vec, propensity_vec_molar
    
    
def eval_dict(to_eval, loc_dict):
    """This function takes a string expression and return the evaluated
    expression using SBML_FUNCT_DICT and the locals() dictionary where
    eval_dict is called.

    Args:
        to_eval (str): the expression to evaluate
        loc_dict (dict): local dictionary from the calling function

    Returns:
        multitype: result of eval command
    """
    return eval(to_eval, loc_dict, SBML_FUNCT_DICT)


def tofloat(val, loc_dict):
    """This function attempts to convert the input val into float

    Args:
        val (str): the expression to evaluate
        loc_dict (dict): local dictionary from the calling function

    Returns:
        float: float equivalent of val
    """
    try:
        return float(val)
    except:
        return float(eval_dict(val, loc_dict))


def is_number(xvar):
    """This function checks if a string xvar is float

    Args:
        xvar (str): input string expression or number

    Returns:
        bool: True if xvar cal be converted to float otherwise False
    """
    try:
        float(xvar)
        return True
    except:
        return False
        
        
def read_topo_data(rfile):
    rows = []
    conc = {}
    with open(rfile, "r") as file:
        last = ""
        for row in file:
            row = row.strip()+" "
            if last == "Function_Definitions":
                if row.strip() != "" and row[0] != "#":
                    exec(row.strip(), locals(), SBML_FUNCT_DICT)
                elif row[0] == "#":
                    last = "#"
            elif last == "#":
                if row.strip() != "" and row[0] != "@":
                    rows.append(row)
                elif row.strip() != "" and row[0] != "#":
                    last = "@"
            elif last == "@":
                if row.strip() != "" and row[0] != "@":
                    cvar = row.split(",")
                    conc[cvar[0].strip()] = tofloat(cvar[1], locals())
            elif row[0] == "#":
                last = "#"
            elif row[0] == "@":
                last = "@"
            elif row.strip().upper() == "FUNCTION_DEFINITIONS:":
                last = "Function_Definitions"
    return rows, conc
    
    
def dict_from_rows(rows):
    ks_dict = {}
    r_dict = {}
    p_dict = {}
    sp_comp = {}
    rxn_rows = len(rows)
    for ih_ind in range(rxn_rows):
        r_dict[ih_ind] = {}
        p_dict[ih_ind] = {}
        col_row = rows[ih_ind].split(":::::")
        row = col_row[0].strip().split(",")

        if len(row) == 3:
            ks_dict[ih_ind] = [
                tofloat(row[1], locals()),
                tofloat(row[2], locals())]
        else:
            ks_dict[ih_ind] = [tofloat(row[1], locals())]
        col_var = row[0].split("<=>")
        if len(col_var) == 1:
            col_var = row[0].split("=>")

        sp_c = col_var[0]
        svar = sp_c.strip().split()
        if len(svar) > 1:
            last = 1
            for xvar in svar:
                if not is_number(xvar) and xvar != "+":
                    r_dict[ih_ind][xvar] = last
                    last = 1
                    if xvar in sp_comp:
                        sp_comp[xvar].add(ih_ind)
                    else:
                        sp_comp[xvar] = {ih_ind}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = svar[0]
            r_dict[ih_ind][xvar] = 1
            if xvar in sp_comp:
                sp_comp[xvar].add(ih_ind)
            else:
                sp_comp[xvar] = {ih_ind}

        sp_c = col_var[1]
        svar = sp_c.strip().split()
        if len(svar) > 1:
            last = 1
            for xvar in svar:
                if (not is_number(xvar) or xvar.lower() == "e") \
                        and xvar != "+":
                    p_dict[ih_ind][xvar] = last
                    last = 1
                    if xvar in sp_comp:
                        sp_comp[xvar].add(ih_ind)
                    else:
                        sp_comp[xvar] = {ih_ind}
                elif is_number(xvar):
                    last = tofloat(xvar, locals())
        else:
            xvar = svar[0]
            p_dict[ih_ind][xvar] = 1
            if xvar in sp_comp:
                sp_comp[xvar].add(ih_ind)
            else:
                sp_comp[xvar] = {ih_ind}

    stoch_var = []

    for sp_c in sp_comp:
        row = []
        for r_ind in range(rxn_rows):
            prod = p_dict[r_ind][sp_c] if sp_c in p_dict[r_ind] else 0
            rect = r_dict[r_ind][sp_c] if sp_c in r_dict[r_ind] else 0
            row.append(prod - rect)
            if len(ks_dict[r_ind]) == 2:
                row.append(-row[-1])
        stoch_var.append(row)
    stoch_var = np.array(stoch_var)
        
    return ks_dict, r_dict, p_dict, sp_comp, stoch_var


class model():
    """This class helps in running BioSANS algorithm as a library"""

    rfile = ""
    miter = 1
    conc_unit = "molar"
    v_volms = 1
    tend = 100
    del_coef = 1.5
    normalize = False
    logx = False
    logy = False
    method = "ODE-2"
    tlen = 1000
    mix_plot = True
    save = False
    out_fname = None
    plot_show = False
    c_input = {}
    vary = ""
    mult_proc = False
    implicit = True
    items = None
    exp_data_file = None
    globals2.init(globals2)
    del_file = None
    del_file2 = None


    def __init__(self, topo=None, sbml=None, FileIn=None, Volume=None):
        """This function initialized the model and model parameters"""

        if not topo:
            if sbml:
                if FileIn.lower() == "molar":
                    sbml_to_topo(sbml, True)
                else:
                    sbml_to_topo(sbml, False)
                topo = sbml+".topo"
            else:
                print("no model file provided")

        try:
            topfile = open(topo, "r")
        except BaseException:
            self.del_file = str(time.time())+"_temp_topo.txt"
            with open(self.del_file, "w") as temp:
                temp.write(topo)
            topo = self.del_file
            topfile = open(topo, "r")

        for row in topfile:
            row = row.strip()+" "
            if row[0] == "#":
                g_g = row.split(",")[1:]
                for xvar in g_g:
                    x_x = [g.strip() for g in xvar.split("=")]
                    if x_x[0] == "Volume":
                        Volume = float(x_x[1])
                    elif x_x[0] == "tend":
                        self.tend = float(x_x[1])
                    elif x_x[0] == "FileUnit":
                        FileIn = x_x[1]
                    elif x_x[0] == "steps":
                        self.tlen = float(x_x[1])

        topfile.close()

        if not FileIn:
            print("concentration unit used in the file not defined : \
                  default to 'molar'")
            FileIn = "molar"
        if not Volume:
            print("Volume not defined : default to 1")
            Volume = 1

        self.rfile = topo
        self.conc_unit = FileIn
        self.v_volms = Volume


    def data(self, exp_data_file=None): 
        """This method instantiate the data file for parameter estimation"""
        if not exp_data_file:
            print("""
                Required experimental data file:

                1) it is a tab delimited file where the first row are
                   headers and succeeding rows are measurements.
                2) the first column is time and succeding columns are
                   components or species
            """)
   
        try:
            expd = open(exp_data_file, "r")
            expd.close()
        except BaseException:
            self.del_file2 = str(time.time())+"_temp_data.txt"
            with open(self.del_file2, "w") as temp:
                glim = "\n"
                exp_datas = exp_data_file.strip().split(glim)
                if len(exp_datas) == 1:
                    glim = "\r\n"
                    exp_datas = exp_data_file.strip().split(glim)             
                exp_data_file = glim.join([xval.strip() for xval in exp_datas])
                temp.write(exp_data_file.replace(",","\t"))
            exp_data_file  = self.del_file2
            
        self.exp_data_file = exp_data_file
        return self


    def plot(self, normalize=False, mix_plot=True, logx=False, logy=False):
        """This function facilitates in the plotting of data"""
        self.normalize = normalize
        self.mix_plot = mix_plot
        self.plot_show = True
        self.logx = logx
        self.logy = logy
        return self


    def save_traj(self, out_fname=None):
        """This function saves the trajectory in to a file"""
        self.save = True
        if not out_fname:
            out_fname = self.rfile+".out.txt"
        self.out_fname = out_fname
        return self


    def extra(self, c_input={}, vary=""):
        """This are extra functions that for now is not needed"""
        self.c_input = c_input
        self.vary = vary
        return self


    def clean(self):
        """This function delete the temporary topology file"""
        if self.del_file:
            delete_file(self.del_file)
            self.del_file = None
        if self.del_file2:
            delete_file(self.del_file2)
            self.del_file2 = None
            
            
    def get_stoich_prop(self, mode="macro"):
        rows, conc = read_topo_data(self.rfile)
        ks_dict, r_dict, p_dict, sp_comp, stoch_var = dict_from_rows(rows)
        if mode == "macro":
            prop = propensity_vec_molar(ks_dict, conc, r_dict, p_dict)
        else:
            prop = propensity_vec(ks_dict, conc, r_dict, p_dict)
        return stoch_var, prop


    def run(self, method=None, ntraj=None, tend=None, step_size_scaler=None,
            steps=None, mult_proc=None, implicit=True, cpu_mult=0.9):
        """This function calls the process function and run the model"""
        globals2.CPU_MULT = cpu_mult
        if method:
            self.method = method
        if ntraj:
            self.miter = ntraj
        if tend:
            self.tend = tend
        if step_size_scaler:
            self.del_coef = step_size_scaler
        if steps:
            self.tlen = steps
        if mult_proc:
            self.mult_proc = mult_proc
        self.implicit = implicit

        data = process(
            rfile=self.rfile,
            miter=self.miter,
            conc_unit=self.conc_unit,
            v_volms=self.v_volms,
            tend=self.tend,
            del_coef=self.del_coef,
            normalize=self.normalize,
            logx=self.logx,
            logy=self.logy,
            method=self.method,
            tlen=self.tlen,
            mix_plot=self.mix_plot,
            save=self.save,
            out_fname=self.out_fname,
            plot_show=self.plot_show,
            c_input=self.c_input,
            vary=self.vary,
            mult_proc=self.mult_proc,
            implicit=self.implicit,
            items=None,
            exp_data_file=self.exp_data_file
        )
        return data