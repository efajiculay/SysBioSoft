"""
                       This is the scrollable_text module

This module controls how text area are placed inside the canvas for
displaying output of computation, equations, etc.

The following is the list of function for this module:

1. save_file
2. tab
3. delete_this
4. prepare_scroll_text

"""

# import sys
# import os
# sys.path.append(os.path.abspath("BioSANS2020"))

from tkinter import Text, INSERT, END, Scrollbar, Frame, Button
from tkinter import filedialog
from sys import platform

from BioSANS2020.myglobal import mglobals as globals2


def save_file(_, text):
    """Opens a file dialog box for saving the constents of the text area"""
    # global file_name
    file = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
    if file is None:
        return -1
    file.write(text.get("0.0", END))
    file.close()
    return "break"


def tab(_, text):
    """insers tab or four spaces in the text area"""
    text.insert(INSERT, " " * 4)
    return "break"


def canvas_update_widgets(_, canvas):
    """This function update the canvas upon adding new objects by moving
    objects to ensure they fit nicely. The arrangement is that the last
    object added appear first (latest goes on top).
    Args:
        canvas : the canvas object"""
    #R = canvas._root().winfo_height()
    #root = canvas._root()
    #width = root.winfo_screenwidth()
    #height = root.winfo_screenheight()
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


def prepare_scroll_text(items):
    """This function prepares the scrollable text area by creating a
    frame as a child of the canvas or items[0]. Scroll capabilities are
    added and a Text area is created as a child of the frame.
    Args:
        items : 3 item list of [canvas, scroll_x, scroll_y]
    Returns:
        text : the text area object where text can be inserted
    """
    # canvas, scroll_x, scroll_y = items
    canvas, _, scroll_y = items
    frame = Frame(canvas, height=455, width=940, borderwidth=10, bd=0)

    hscroll = Scrollbar(frame, orient="horizontal")
    hscroll.pack(side="bottom", fill="x")
    vscroll = Scrollbar(frame)
    vscroll.pack(side="right", fill="y")

    font_val = 11 if platform.lower() != "darwin" else 15
    text = Text(frame, width=int(106 * canvas.winfo_width() / 972), height=25,
                fg='blue', font=("Courier New", font_val, "bold"), wrap="none")
    text.configure(xscrollcommand=hscroll.set, yscrollcommand=vscroll.set)
    text.pack(side="top", fill="both", expand=True)

    frame.pack(side="top", fill="x", expand=True)

    vscroll.config(command=text.yview)
    hscroll.config(command=text.xview)

    fframe = canvas.create_window(
        0, 450 * globals2.PLOT_I, anchor='nw', window=frame)
    bttn = Button(frame, text=" X ", fg='red', highlightcolor='blue',
                  bg='white', height=1, relief='raised',
                  command=lambda: delete_this(fframe, canvas))
    bttn.place(rely=0.0, relx=1.0, x=-15, y=0, anchor="ne")

    # canvas.configure(yscrollcommand=scroll_y.set,
    #                  xscrollcommand=scroll_x.set)
    canvas.configure(yscrollcommand=scroll_y.set)
    canvas.configure(scrollregion=canvas.bbox("all"))
    globals2.PLOT_I = globals2.PLOT_I + 1
    canvas_update_widgets(None, canvas)

    text.bind("<Button-2>", lambda e: save_file(e, text))
    text.bind("<Tab>", lambda e: tab(e, text))
    # text.bind("<Button-3>",lambda e: print("kkkk"))
    canvas.bind("<Configure>", lambda e: canvas_update_widgets(e, canvas))
    # canvas.bind("<Configure>", lambda e: [ canvas.itemconfig(xx,width= \
    #     canvas.winfo_width()) for xx in canvas.find_all() ])
    return text
