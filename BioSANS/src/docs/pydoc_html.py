import os

summary = {}
for dirpath, dirs, files in os.walk('BioSANS2020'): 
	for filename in files:
		fname = os.path.join(dirpath,filename)
		if fname.endswith('.py'):
			name = dirpath.replace("\\",".")+"."+filename.replace(".py","")
			os.system("pydoc -w "+name)
			

   