from BioSANS2020.prepcodes.process import process
from BioSANS2020.model.fileconvert.process_sbml import process_sbml as sbml_to_topo
from BioSANS2020.myglobal import mglobals as globals2
#from BioSANS2020.myglobal import proc_global as proc_global


class model():

	globals2.init(globals2)
	#proc_global.init(proc_global)

	def __init__(
		self,
		topo=None,
		sbml=None,
		ntraj=1,
		FileIn=None,
		Volume=None,
		tend=None,
		steps=None,
		step_size_scaler=10,
		normalize=False,
		mix_plot= True,
		save=False,
		out_fname=None,
		plot_show=False,
		c_input={},
		vary="",
		mult_proc=False,
		cpu_mult=0.9,
		implicit=True,
		exp_data_file=None
	):
	
		if not topo:
			if FileIn.lower() == "molar":			
				sbml_to_topo(sbml,True)
			else:
				sbml_to_topo(sbml,False)
			topo = sbml+".topo"
	
		topfile = open(topo, "r")
		for row in topfile:
			if row[0] == "#":
				g_g = row.split(",")[1:]
				for xvar in g_g:
					x_x = [g.strip() for g in xvar.split("=")]
					if x_x[0] == "Volume":
						Volume = x_x[1]
					elif x_x[0] == "tend":
						tend = x_x[1]
					elif x_x[0] == "FileUnit":
						FileIn = x_x[1]
					elif x_x[0] == "steps":
						steps = x_x[1]
			
		if not FileIn:
			print("concentration unit used in the file not defined : default = 'molar'")
			FileIn = "molar"
		if not Volume:
			print("Volume not defined : default = 1")
			Volume = 1
		if not tend:
			tend = 100
			print("end time of simulation not defined : default = 100")
		if not steps:
			print("number of steps to take not defined : default = 1000")
			steps = 1000
		if not out_fname:
			out_fname = topo+".out.txt"
		
	
		self.rfile    	= topo
		self.miter      = ntraj
		self.conc_unit	= FileIn
		self.v_volms 	= Volume
		self.tend       = tend 
		self.del_coef	= step_size_scaler
		self.normalize	= normalize
		self.logx       = False
		self.logy       = False
		self.tlen       = steps
		self.mix_plot	= True
		self.save       = save
		self.out_fname	= out_fname
		self.plot_show	= False
		self.c_input    = {}
		self.vary       = ""
		self.mult_proc	= mult_proc
		self.implicit   = implicit
		self.exp_data_file = exp_data_file
		globals2.CPU_MULT = cpu_mult
	
	
	def run(self, method):
		return process(
			rfile=self.rfile,
			miter=self.miter,
			conc_unit=self.conc_unit,
			v_volms=self.v_volms,
			tend=self.tend,
			del_coef=self.del_coef,
			normalize=self.normalize,
			logx=self.logx,
			logy=self.logy,
			method=method,
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