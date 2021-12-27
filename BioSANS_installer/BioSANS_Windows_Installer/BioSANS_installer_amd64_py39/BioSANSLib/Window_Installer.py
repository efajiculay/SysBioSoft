import sys
import os
import subprocess
from subprocess import Popen
import pkg_resources

#print(os.environ["HOMEPATH"])
#print(sys.version)

python_ok = False
if sys.version_info[0] == 3 and sys.version_info[1] >= 7:
	python_ok = True

pythonCommand = "C:/BioSANS2021/python39amd64/Python39/python"
pipCommand = pythonCommand+" -m pip"

if python_ok:

	try:
		os.system(pipCommand+" "+"install --upgrade pip")
	except:
		pass
	
	com = pipCommand+" "+"install BioSANS2020 --upgrade"
	rs = Popen(com).communicate()[0]
	
	Drive = os.path.abspath(os.environ["HOMEPATH"]).split(":")[0]+":"
	try:
		os.mkdir(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS").replace(Drive,"C:"),0o777)
	except:
		os.remove(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS/BioSANS.bat"))
		os.remove(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS/BioSANS.ico"))
		
	Shcut = open(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS/BioSANS.bat").replace(Drive,"C:"),"w")
	Shcut.write("@echo off\n")
	Shcut.write("start "+pythonCommand+"w -m BioSANS2020.BioSANS")
	Shcut.close()
	
		
else:
	print("Please install python 3.7 or higher first")
	
