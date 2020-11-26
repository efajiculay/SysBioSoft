from tkinter import Canvas, Scrollbar
#frame1a = Frame(frame1, height = 510, width = 1045, bg='#8c8c8c', borderwidth=2)
#frame1a.place(x=2,y=75) 

def prepare_frame_for_plot(frame1a, width=1030, height=465):
	canvas = Canvas(frame1a, width=width, height=height)
	canvas.grid(row=0, column=0)	
	scroll_x = Scrollbar(frame1a, orient="horizontal", command=canvas.xview)
	scroll_x.grid(row=1, column=0, sticky="ew")
	scroll_y = Scrollbar(frame1a, orient="vertical", command=canvas.yview)
	scroll_y.grid(row=0, column=1, sticky="ns")
	canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
	canvas.configure(scrollregion=canvas.bbox("all"))	
	return (canvas, scroll_x, scroll_y)
