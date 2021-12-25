"""

                   This module is the mglobals module

Variables global to  all modules are listed here and instantiated in one
of the main calling module such as BioSANS module and BioSSL module.

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

PLOT_I = 0
INT_VARS = []
CONTAINER = []
PLOTTED = []
MODIFIED = {}
PROP_MODIFIED = {}
CON_BOUNDARY = {}
TO_CONVERT = ""
SETTINGS = {}
EXEC_FUNCTIONS = []
DELAY_LIST = {}
TCHECK = []


def init(self):
    """Re initiate global variable just in case needed
    """
    self.PLOT_I = 0
    self.INT_VARS = []
    self.CONTAINER = []
    self.PLOTTED = []
    self.MODIFIED = {}
    self.PROP_MODIFIED = {}
    self.CON_BOUNDARY = {}
    self.TO_CONVERT = ""
    self.SETTINGS = {}
    self.EXEC_FUNCTIONS = []
    self.DELAY_LIST = {}
    self.TCHECK = []
