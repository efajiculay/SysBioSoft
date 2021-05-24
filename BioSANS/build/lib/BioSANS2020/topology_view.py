import sys
import os
sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.scrollable_text import *

def view_topo(topo,items):	

	text = prepare_scroll_text(items)
	ffprint = lambda x: text.insert(INSERT," ".join([str(y) for y in x]))	
	with open(topo) as f:
		for x in f:
			ffprint([x.replace("\t"," "*4)])
	return text