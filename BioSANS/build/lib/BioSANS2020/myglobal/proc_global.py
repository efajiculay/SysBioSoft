"""

                 This module is the proc_global module

The  main  purpose of this module is  to observed module wide processes.

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

import multiprocessing as mp


MANAGER = None
LST = None


def init(self):
    """Initiate process related global variables
    """
    self.MANAGER = mp.Manager()
    self.LST = self.MANAGER.list()
