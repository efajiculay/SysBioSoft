#import sys
#import os
# sys.path.append(os.path.abspath("BioSANS2020"))

def init():
    global process_call, plot_i, IntVars, Container, plotted, modified, Constant, PropModified, conBoundary, settings, execFunctions, delay_list
    process_call = 0
    plot_i = 0
    IntVars = []
    Container = []
    plotted = []
    modified = {}
    Constant = {}
    PropModified = {}
    conBoundary = {}
    toConvert = ""
    settings = {}
    execFunctions = []
    delay_list = {}
