"""

				This is the topology_view module

	This module helps to load topology files into a text area for
displaying in BioSANS.

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text
from BioSANS2020.gui_functs.scrollable_text import INSERT


def view_topo(topo, items):
    """This function open topology files and display the topology into a
    text area for further simulation or for editing.

    Args:
        topo (string): topology file name
        items (tuple): (canvas, scroll_x, scroll_y)

    Returns:
        tkinter.Text: text area where topology is displayed. A user can
            type and modify the topology
    """
    text = prepare_scroll_text(items)

    def ffprint(xvar):
        return text.insert(INSERT, " ".join([str(y) for y in xvar]))
    with open(topo) as ffile:
        for xvar in ffile:
            ffprint([xvar.replace("\t", " " * 4)])
    return text
