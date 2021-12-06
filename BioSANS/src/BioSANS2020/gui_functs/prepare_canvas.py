"""

           This is the prepare_canvas module

This module prepare the canvas for drawing figures, plots and images.

The following is the list of function for this module:

1. prepare_frame_for_plot

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from tkinter import Canvas, Scrollbar
# frame1a = Frame(frame1, height = 510, width = 1045, bg='#8c8c8c',
#                 borderwidth=2)
# frame1a.place(x=2,y=75)


def prepare_frame_for_plot(frame1a, width=1030, height=465):
    """This function prepares the frame that will serve as a parent of
    the canvas that will serve as a parent of many frames that will
    contain figures, image, and plots.
    Args:
        frame1a : the parent frame
        width : width of frame1a
        height : height frame1a
    Returns:
        tuple : (canvas, scroll_x, scroll_y)

    """
    canvas = Canvas(frame1a, width=width, height=height)
    canvas.pack(side="left", fill="both", expand=True)
    scroll_x = Scrollbar(frame1a, orient="horizontal", command=canvas.xview)
    # scroll_x.pack(side = "bottom", fill = "x")
    scroll_y = Scrollbar(frame1a, orient="vertical", command=canvas.yview)
    scroll_y.pack(side="right", fill="y")
    # canvas.configure(yscrollcommand=scroll_y.set,
    #                  xscrollcommand=scroll_x.set)
    canvas.configure(yscrollcommand=scroll_y.set)
    # canvas.configure(scrollregion=canvas.bbox("all"))
    return (canvas, scroll_x, scroll_y)
