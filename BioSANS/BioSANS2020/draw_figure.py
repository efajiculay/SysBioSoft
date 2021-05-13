from tkinter import Frame, Canvas, Checkbutton, TOP, BOTH, IntVar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk as NavigationToolbar2TkAgg
import mglobals as globals2

def draw_figure(items,figure,loc=(0, 0)):
	if items:
		canvas,scroll_x,scroll_y = items
		canva = Canvas(canvas,height = 426, width = 1030,bg='#ccffcc')	
		frame = Frame(canva, height = 425, width = 1000, borderwidth=0,bd=0)
		frame.pack(side='top')		
		figure_canvas = FigureCanvasTkAgg(figure, frame)
		canva.pack()
		figure_canvas.get_tk_widget().pack(side='top')
		toolbar = NavigationToolbar2TkAgg(figure_canvas, frame)
		toolbar.update()
		figure_canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)	
		globals2.IntVars.append(IntVar(value=-1))
		C1 = Checkbutton(canva, text = "", variable = globals2.IntVars[-1], onvalue = globals2.plot_i, offvalue = -1, height=2, width = 2,bg="#ccffcc")
		C1.place(x=1000,y=200)		
		wind1 = canva.create_window(2, 2, anchor='nw', window=frame)
		wind2 = canva.create_window(1020, 216, anchor='ne', window=C1,tags="gg")
		canvas.create_window(0, 426*globals2.plot_i, anchor='nw', window=canva)
		globals2.Container.append([canvas,wind1])
		#canvas.move(wind1,100,100)
		canvas.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
		canvas.configure(scrollregion=canvas.bbox("all"))	
		globals2.plot_i = globals2.plot_i+1	
	
	
	