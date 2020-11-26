import multiprocessing as mp

def init():
	global manager, lst
	manager = mp.Manager()	
	lst = manager.list()
	