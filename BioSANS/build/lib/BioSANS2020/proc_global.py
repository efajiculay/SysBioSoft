import sys
import os
sys.path.append(os.path.abspath("BioSANS2020"))

import multiprocessing as mp

def init():
	global manager, lst
	manager = mp.Manager()	
	lst = manager.list()
	