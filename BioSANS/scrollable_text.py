from tkinter import Text, INSERT, END, Scrollbar, RIGHT, LEFT, Y, Canvas, VERTICAL, Frame, Menu
import mglobals as globals2
from tkinter import filedialog

def save_file(e,text):
	global file_name
	file = filedialog.asksaveasfile(mode='w', defaultextension=".txt")
	if file is None:
		return
	file.write(text.get("0.0",END))
	file.close()
	return "break"
	
def prepare_scroll_text(items):
	canvas,scroll_x,scroll_y = items
	frame = Frame(canvas, height = 415, width = 940, borderwidth=0,bd=0)
	frame.pack(side='top')	
	text = Text(frame,width=106,height=24,fg='blue',font = ("Courier New", 11, "bold"))
	text.pack(side=LEFT)
	tscroll = Scrollbar(frame,orient=VERTICAL)
	tscroll.pack(side=RIGHT, fill=Y)
	text.configure(yscrollcommand=tscroll.set)
	tscroll.config(command=text.yview)
	canvas.create_window(0, 426*globals2.plot_i, anchor='nw', window=frame)
	canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
	canvas.configure(scrollregion=canvas.bbox("all"))
	globals2.plot_i = globals2.plot_i + 1
	
	text.bind("<Button-2>",lambda e: save_file(e,text) )
	return text
