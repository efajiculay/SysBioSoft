#import sys
#import os
#sys.path.append(os.path.abspath("BioSANS2020"))

from tkinter import Text, INSERT, END, Scrollbar, RIGHT, LEFT, Canvas, Frame, Menu, Widget
from BioSANS2020.myglobal import mglobals as globals2
from tkinter import filedialog
from sys import platform

def save_file(e,text):
	global file_name
	file = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
	if file is None:
		return
	file.write(text.get("0.0",END))
	file.close()
	return "break"
	
def tab(e,text):
    text.insert(INSERT, " " * 4)
    return "break"

#https://www.it.uu.se/edu/course/homepage/scriptprog/st07/material/tkinterref	
def canvas_update_widgets(e,canvas):
	#R = canvas._root().winfo_height()
	#root = canvas._root()
	#width = root.winfo_screenwidth()
	#height = root.winfo_screenheight()	
	#print(canvas.winfo_children())
	H = canvas.winfo_height()
	W = canvas.winfo_width()
	objCan = canvas.find_all()
	for x in objCan:
		canvas.itemconfig(x,height=H,width=W-5)
		x1, y1, x2, y2 = canvas.bbox(x)
		canvas.move(x,0,(x-1)*H-y1+3)
	canvas.configure(scrollregion=canvas.bbox("all"))
	return "break"

def prepare_scroll_text(items):
	canvas,scroll_x,scroll_y = items
	frame = Frame(canvas, height = 455, width = 940, borderwidth=10,bd=0)
	
	Hscroll = Scrollbar(frame, orient = "horizontal")	
	Hscroll.pack(side = "bottom", fill = "x")
	Vscroll = Scrollbar(frame)
	Vscroll.pack(side = "right", fill = "y")
	
	font_val = 11 if platform.lower() != "darwin" else 15
	text = Text(frame,width=int(106*canvas.winfo_width()/972),height=25,fg='blue',font = ("Courier New", font_val, "bold"),wrap="none")
	text.configure(xscrollcommand=Hscroll.set,yscrollcommand=Vscroll.set)
	text.pack(side="top",fill="both",expand=True)
	
	frame.pack(side="top",fill="x",expand=True)	

	Vscroll.config(command=text.yview)
	Hscroll.config(command=text.xview)
	
	canvas.create_window(0, 450*globals2.plot_i, anchor='nw', window=frame)
	#canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
	canvas.configure(yscrollcommand=scroll_y.set)
	canvas.configure(scrollregion=canvas.bbox("all"))
	globals2.plot_i = globals2.plot_i + 1
	canvas_update_widgets(None,canvas)
	
	text.bind("<Button-2>",lambda e: save_file(e,text) )
	text.bind("<Tab>",lambda e: tab(e,text))
	canvas.bind("<Configure>", lambda e: canvas_update_widgets(e,canvas))
	#canvas.bind("<Configure>", lambda e: [ canvas.itemconfig(xx,width=canvas.winfo_width()) for xx in canvas.find_all() ])
	return text