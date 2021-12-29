import time
from os import remove as delete_file
from BioSANS2020.prepcodes.process import process
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2



class model():

	rfile=""
	miter=1
	conc_unit="molar"
	v_volms=1
	tend=100
	del_coef=1.5
	normalize=False
	logx=False
	logy=False
	method="ODE-2"
	tlen=1000
	mix_plot=True
	save=False
	out_fname=None
	plot_show=False
	c_input={}
	vary=""
	mult_proc=False
	implicit=True
	items=None
	exp_data_file=None
	globals2.init(globals2)
	del_file = None

	def __init__(
		self,
		topo=None,
		sbml=None,
		FileIn=None,
		Volume=None
	):
	
		if not topo:
			if sbml:
				if FileIn.lower() == "molar":			
					sbml_to_topo(sbml,True)
				else:
					sbml_to_topo(sbml,False)
				topo = sbml+".topo"
			else:
				print("no model file provided")
	
		try:
			topfile = open(topo, "r")
		except:
			self.del_file = str(time.time())+"_temp_topo.txt"
			with open(self.del_file, "w") as temp:
				temp.write(topo)
			topo = self.del_file
			topfile = open(topo, "r")
		
		for row in topfile:
			row = row.strip()+" "
			if row[0] == "#":
				g_g = row.split(",")[1:]
				for xvar in g_g:
					x_x = [g.strip() for g in xvar.split("=")]
					if x_x[0] == "Volume":
						Volume = float(x_x[1])
					elif x_x[0] == "tend":
						self.tend = float(x_x[1])
					elif x_x[0] == "FileUnit":
						FileIn = x_x[1]
					elif x_x[0] == "steps":
						self.tlen = float(x_x[1])
						
		topfile.close()
			
		if not FileIn:
			print("concentration unit used in the file not defined : default to 'molar'")
			FileIn = "molar"
		if not Volume:
			print("Volume not defined : default to 1")
			Volume = 1
	
		self.rfile    	= topo
		self.conc_unit	= FileIn
		self.v_volms 	= Volume 
		
		
	def data(self, exp_data_file=None):
		if not exp_data_file:
			print("""
				Required experimental data file:
				
				1) it is a tab delimited file where the first row are
				   headers and succeeding rows are measurements.
				2) the first column is time and succeding columns are
				   components or species
			""")
		self.exp_data_file = exp_data_file
		return self
		
		
	def plot(self, normalize=False, mix_plot=True, logx=False, logy=False):
		self.normalize	= normalize
		self.mix_plot	= mix_plot
		self.plot_show	= True
		self.logx       = logx
		self.logy       = logy
		return self
		
		
	def save_traj(self, out_fname=None):
		self.save = True
		if not out_fname:
			out_fname = self.rfile+".out.txt"
		self.out_fname	= out_fname
		return self
		
		
	def extra(self, c_input={}, vary=""):
		self.c_input = c_input
		self.vary    = vary
		return self
		
		
	def clean(self):
		if self.del_file:
			delete_file(self.del_file)
			self.del_file = None
			
	
	def run(self, method=None, ntraj=None, tend=None, step_size_scaler=None,
	        steps=None, mult_proc=None, implicit=True, cpu_mult=0.9):
		globals2.CPU_MULT = cpu_mult
		if method:
			self.method = method
		if ntraj:
			self.miter = ntraj
		if tend:
			self.tend = tend 
		if step_size_scaler:
			self.del_coef = step_size_scaler
		if steps:
			self.tlen = steps
		if mult_proc:
			self.mult_proc = mult_proc
		self.implicit = implicit

		data = process(
			rfile=self.rfile,
			miter=self.miter,
			conc_unit=self.conc_unit,
			v_volms=self.v_volms,
			tend=self.tend,
			del_coef=self.del_coef,
			normalize=self.normalize,
			logx=self.logx,
			logy=self.logy,
			method=self.method,
			tlen=self.tlen,
			mix_plot=self.mix_plot,
			save=self.save,
			out_fname=self.out_fname,
			plot_show=self.plot_show,
			c_input=self.c_input,
			vary=self.vary,
			mult_proc=self.mult_proc,
			implicit=self.implicit,
			items=None,
			exp_data_file=self.exp_data_file
		)
		return data