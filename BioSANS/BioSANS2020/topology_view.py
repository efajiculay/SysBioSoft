from scrollable_text import *

def view_topo(topo,items):	

	text = prepare_scroll_text(items)
	ffprint = lambda x: text.insert(INSERT," ".join([str(y) for y in x]))	
	with open(topo) as f:
		for x in f:
			ffprint([x])
	return text