"""
                       This is the draw_figure module

This module controls how plots are drawn inside BioSANS canvas

The following is the list of function for this module:

1. canvas_update_widgets
2. delete_this
3. draw_figure

"""


# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

# from tkinter import Frame, Canvas, Checkbutton, IntVar, Button
from tkinter import Frame, Canvas, IntVar, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, \
    NavigationToolbar2Tk as NavigationToolbar2TkAgg
from BioSANS2020.myglobal import mglobals as globals2


def canvas_update_widgets(_, canvas):
    """This function update the canvas upon adding new objects by moving
    objects to ensure they fit nicely. The arrangement is that the last
    object added appear first (latest goes on top).
    Args:
        canvas : the canvas object"""
    # R = canvas._root().winfo_height()
    # root = canvas._root()
    # width = root.winfo_screenwidth()
    # height = root.winfo_screenheight()
    # print(canvas.winfo_children())
    canvas_height = canvas.winfo_height()
    canvas_width = canvas.winfo_width()

    obj_canv = canvas.find_all()
    obj_len = len(obj_canv)
    ind = 0
    for xvar in obj_canv:
        canvas.itemconfig(xvar, height=canvas_height, width=canvas_width - 5)
        # xx1, yy1, xx2, yy2 = canvas.bbox(xvar)
        _, yy1, _, _ = canvas.bbox(xvar)
        canvas.move(xvar, 0, (obj_len - ind) * canvas_height - yy1 + 3)
        ind = ind + 1
    canvas.configure(scrollregion=canvas.bbox("all"))
    return "break"


def delete_this(frame, canvas):
    """This function removes objects from the canvas and reorder objects"""
    canvas.delete(frame)
    canvas_update_widgets(None, canvas)


# def draw_figure(items, figure, loc=(0, 0)):
def draw_figure(items, figure, _):
    """This function draws figure to the canvas by first adding a canvas
    object to the canvas, adding a frame to the newly created canvas and
    putting the figure to FigureCanvasTkAgg with frame as the parent.
    Args:
        items : three dimensional array of [canvas, scroll_x, scroll_y]
        figure: the figure to be drawn in canvas i.e. plt.gcf() figure
    """
    if items:
        canvas, scroll_x, scroll_y = items
        canva = Canvas(canvas, height=426, width=1030, bg='#ccffcc')
        frame = Frame(canva, height=425, width=1000, borderwidth=0, bd=0)
        frame.pack(side='top')
        figure_canvas = FigureCanvasTkAgg(figure, frame)
        canva.pack(fill="both", expand=True)
        figure_canvas.get_tk_widget().pack(side='top')
        toolbar = NavigationToolbar2TkAgg(figure_canvas, frame)
        toolbar.update()
        figure_canvas._tkcanvas.pack(side="top")
        globals2.IntVars.append(IntVar(value=-1))
        # C1 = Checkbutton(canva, text = "", variable = globals2.IntVars[-1],
        #                  onvalue = globals2.plot_i, offvalue = -1,
        #                  height=2, width = 2,bg="#ccffcc")
        # C1.place(x=1000,y=200)
        wind1 = canva.create_window(2, 2, anchor='nw', window=frame)
        # wind2 = canva.create_window(1020, 216, anchor='ne', window=C1,
        #                             tags="gg")
        fframe = canvas.create_window(
            0, 450 * globals2.plot_i, anchor='nw', window=canva)
        buttn = Button(frame, text=" X ", fg='red', highlightcolor='blue',
                       bg='white', height=1, relief='raised',
                       command=lambda: delete_this(fframe, canvas))
        buttn.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")

        globals2.Container.append([canvas, wind1])
        # canvas.move(wind1,100,100)
        canvas.configure(yscrollcommand=scroll_y.set,
                         xscrollcommand=scroll_x.set)
        canvas.configure(scrollregion=canvas.bbox("all"))
        globals2.plot_i = globals2.plot_i + 1
        canvas_update_widgets(None, canvas)
