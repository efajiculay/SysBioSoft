#import sys
#import os
#sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import *

def new_file(items):	
	text = prepare_scroll_text(items)
	return text