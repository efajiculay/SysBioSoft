"""

                    This module is the new_file module

    This module helps in displaying a new text area for typing topology,
codes, or regular text.

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from BioSANS2020.gui_functs.scrollable_text import prepare_scroll_text


def new_file(items):
    """This function returns a blank text area for typing
    Args:
        items (tuple): (canvas, scroll_x, scroll_y)

    Returns:
        tkinter.Text: text area where a user can type
    """
    text = prepare_scroll_text(items)
    return text
