import sys
import os
import pathlib
from process import *
from ssl_calls import *
import pandas as pd
import mglobals as globals
import matplotlib.pyplot as plt
globals.init()

cwd = str(pathlib.Path().absolute())
trj = {}
try:
	os.mkdir("SSLtemps")
except:
	pass

def get_input():
	try:
		return raw_input("> ")
	except:
		return input("> ")
		
#propagate 
#A => B, 0.2 & B => C, 0.3 
#where A=100 & B=0.2 & C=0 
#using CLE 
#with tn=50 & tlen=1000 & miter=2 & mult_proc=True & fout=Traj1 & plot=False & mixp=False;

#propagate A => B, 0.2 & B => C, 0.3 where A=100 & B=0.2 & C=0 using LNA with tn=50 & tlen=1000;
#propagate A => B, 0.2 & B => C, 0.3 where A=100 & B=0.2 & C=0 using LNA2 with tn=50 & tlen=1000;	

def process_command(command):
	global cwd, trj
	rowc = command.strip().split()
	miter = 1
	tn = 100
	Vol = 1.0e-10
	tsc = 10
	tlen=100
	mult_proc = False
	fout = "temp_traj"
	plot=True
	mixp=True
	if rowc[0].lower() == "propagate":
		rlen = len(rowc)
		i = 1
		rxns = ""
		while rowc[i] != "where" and i<rlen:
			rxns = rxns+rowc[i]
			i = i+1
		rxns = rxns.split("&")
		
		fFile = cwd+"/SSLtemps/temp.txt"
		f = open(fFile,"w")
		f.write("#REACTIONS\n")
		for rx in rxns:
			f.write(rx+"\n")
		f.write("\n")
		
		Mname = ""
		f.write("@CONCENTRATION\n")
		conc = ""
		
		i = i+1
		while rowc[i] != "using" and i<rlen:
			conc = conc+rowc[i]
			Mname = rowc[i+2]
			i = i+1

		for cx in conc.split("&"):
			f.write(cx.replace("=",",")+"\n")
		f.write("\n")			
		f.close()
		
		i = i+3
		opts = ""
		while i<rlen:
			opts = opts+rowc[i]
			i = i+1
		opts = opts.split("&")
		
		optsv = {}
		for op in opts:
			opr = op.split("=")
			optsv[opr[0].strip()] = opr[1]
			
		if 'tn' in optsv:
			tn = float(optsv['tn'])
		if 'tlen' in optsv:
			tlen = int(optsv['tlen'])
		if 'Vol' in optsv:
			Vol = float(optsv['Vol'])
		if 'tsc' in optsv:
			tsc = float(optsv['tsc'])				
		if 'miter' in optsv:
			miter = int(optsv['miter'])	
		if 'mult_proc' in optsv:
			mult_proc = optsv['mult_proc'].lower() == "true"
		if 'fout' in optsv:
			fout = optsv['fout']
		if 'plot' in optsv:
			plot = optsv['plot'].lower() == "true"		
		if 'mixp' in optsv:
			mixp = optsv['mixp'].lower() == "true"				
			
		Fname = cwd+"/SSLtemps/"+fout
		process(
			rfile=fFile,
			miter=miter,
			inMolar=False,
			Vm=Vol,
			tn=tn,
			delX=tsc,
			normalize=False,
			logx=False,
			logy=False,
			method=Mname,
			tlen=tlen,
			mix_plot=mixp,
			save=True,
			out_fname=Fname,
			plot_show=plot,
			Cinput={},
			vary = "",
			mult_proc=mult_proc,
			items=None
		)
	elif rowc[0].lower() == "load":
		try:
			sslfile = cwd+"/SSLtemps/"+rowc[1]
			f = open(sslfile)	
		except:
			sslfile = rowc[1]
			f = open(sslfile)
		command = ""
		for row in f:
			if len(row.strip())>0:
				if row.strip()[-1] == ";":
					command = command + row.strip().replace(";","")
					process_command(command)
					command = ""
				else:
					command = command + row.strip() + " " 
		f.close()	
	elif rowc[0].lower() == "pwd":
		print(cwd)
	elif rowc[0].lower() == "mkdir":
		os.mkdir(rowc[1])
	elif rowc[0].lower() == "ls":
		dirs = []
		try:
			if len(rowc)>1:
				dirs = os.listdir(cwd+"/"+rowc[1])
			else:
				dirs = os.listdir(cwd)
		except:
			pass
		for x in dirs:
			if os.path.isdir(x):
				print("directory : "+x)
			else:
				print("file      : "+x)
	elif rowc[0].lower() == "cd":
		try:
			abspath = os.path.abspath(cwd+"/"+rowc[1])
			dname = os.path.dirname(abspath)
			os.chdir(dname)
			cwd = str(abspath)
		except:
			print("cannot change dir")
	elif rowc[0].lower() == "read_traj":
		if len(rowc) == 4:
			if rowc[2].lower() == "as":
				name = rowc[3].strip()
				try:
					sslfile = cwd+"/SSLtemps/"+rowc[1]
					trj[name] = load_data_traj(sslfile)	
				except:
					sslfile = rowc[1]	
					trj[name] = load_data_traj(sslfile)				
			else:
				print("use as to assign content to variable")
				print("read_traj filename as variable")
		else:
			print("Invalid syntax")
	elif rowc[0].lower() == "print":
		try:
			if len(rowc) == 3:
				try:
					print(trj[rowc[1].strip()][eval(rowc[2].strip())])
				except:
					print(trj[rowc[1].strip()][eval("'"+rowc[2].strip()+"'")])
			else:
				print(trj[rowc[1].strip()])
		except:
			pass
	elif rowc[0].lower() == "plot":
		if len(rowc) == 4:
			trj[rowc[1].strip()].plot(x=rowc[2],y=rowc[3],kind='scatter',s=1)
			plt.show()
		elif len(rowc) > 4:
			qx = rowc[3:]
			plt.xlabel(rowc[2])
			for x in qx:
				plt.scatter(trj[rowc[1]][rowc[2]],trj[rowc[1]][x],label=x,s=1)
			plt.legend()
			plt.show()			
		else:
			print("Invalid syntax")
	elif rowc[0].lower() == "calc_covariance":
		try:
			calc_covariance(trj[rowc[1].strip()], int(rowc[2]))
		except:
			calc_covariance(trj[rowc[1].strip()], 100)
	elif rowc[0].lower() == "prob_density":
		Prob_density(trj[rowc[1].strip()],"SSLtemps/"+rowc[1].strip())
	elif rowc[0].lower() == "prob_density_wtime":
		Prob_density_wtime(trj[rowc[1].strip()],"SSLtemps/"+rowc[1].strip(),"Mname")
	elif rowc[0].lower() == "calc_average":
		try:
			calc_average_conc_at_tend(trj[rowc[1].strip()], int(rowc[2]))
		except:
			calc_average_conc_at_tend(trj[rowc[1].strip()], 100)
	elif rowc[0].lower() == "length":
		try:
			print(len(trj[rowc[1].strip()]))
		except:
			pass
	elif rowc[0].lower() == "pdread_traj":
		if len(rowc) == 4:
			if rowc[2].lower() == "as":
				name = rowc[3].strip()
				try:
					sslfile = cwd+"/SSLtemps/"+rowc[1]
					trj[name] = pd.read_csv(sslfile,delimiter="\t")	
				except:
					sslfile = rowc[1]	
					trj[name] = pd.read_csv(sslfile,delimiter="\t")		
			else:
				print("use as to assign content to variable")
				print("read_traj filename as variable")
		else:
			print("invalid syntax")
if __name__ == '__main__':
	print("###############################################################")
	print("Welcome to Biostoch-SSL commandline interface\n") 
		
	row = " "
	command = " "

	while row.strip() != "quit":
		row = " "+get_input()
		if row[-1]==";":
			command = command + " " + row.strip().replace(";","")
			if command.strip() != "":
				process_command(command)
			command = " "
			row = " "
		else:
			command = command + " " + row




