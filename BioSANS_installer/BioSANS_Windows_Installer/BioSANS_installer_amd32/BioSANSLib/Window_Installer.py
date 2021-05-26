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

pythonCommand = "C:/BioSANS2021b/Python39/python"
pipCommand = pythonCommand+" -m pip"

if python_ok:

	try:
		os.system(pipCommand+" "+"install --upgrade pip")
	except:
		pass

	"""
	required = {
		"DateTime",
		"pathlib",
		"Pillow",
		"numpy",
		"thread6",
		"pandas",
		"matplotlib==3.3.3",
		"python-libsbml",
		"scipy",
		"sdeint",
		"func_timeout",
		"sympy"
	}
	
	installed = {pkg.key for pkg in pkg_resources.working_set}
	missing = required - installed
	
	for item in missing:	
		try:
			os.system(pipCommand+" "+"install "+item+"")
		except:
	""" 
														
		
	com = pipCommand+" "+"install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple BioSANS2020"
	rs = Popen(com).communicate()[0]
	
	"""
	A = subprocess.check_output([pipCommand, 'show', 'BioSANS2020'])
	A = str(A).split("\\n")
	install_dir = ""
	for row in A:
		line = row.split(":")
		if line[0].strip() == "Location":
			install_dir = ("".join(line[1:]).strip("\\r\\n").replace("c\\","c:/").replace("\\","/").replace("//","/"))
		
	if install_dir != "":
		install_dir = str(install_dir)
		print()
		print("To run BioSANS GUI, type the following commands")
		print()
		print(pythonCommand+" "+install_dir+"/BioSANS2020/BioSANS.py")
		print()
		print()
		print("To run BioSANS SSL console, type the following commands")
		print()
		print(pythonCommand+" "+install_dir+"/BioSANS2020/BioSSL.py")
		print()

		Drive = os.path.abspath(os.environ["HOMEPATH"]).split(":")[0]+":"
		try:
			os.mkdir(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS").replace(Drive,"C:"),0o777)
		except:
			pass
		Shcut = open(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS/BioSANSb.bat").replace(Drive,"C:"),"w")
		Shcut.write("@echo off\n")
		Shcut.write("start "+pythonCommand+"w "+install_dir+"/BioSANS2020/BioSANS.py")
		Shcut.close()
		print()
		print("You may also launched BioSANS by double clicking the BioSANS.bat file in the desktop")
		print()	
		
		p = Popen("short_cut_create.bat", cwd=os.getcwd())
		stdout, stderr = p.communicate()
		
	else:
		print("software installed but directory of installation not automatically grabbed")
		print("Locate installation directory using\n")
		print(pipCommand+" "+"show BioSANS2020")
		print()
		print()
		print("To run BioSANS GUI, type the following commands")
		print()
		print(pythonCommand+" "+install_dir+"/BioSANS2020/BioSANS")
		print()
		print("To run BioSANS SSL console, type the following commands")
		print()
		print(pythonCommand+" "+install_dir+"/BioSANS2020/BioSSL.py")
		print()
	"""
	
	print("####################################################################")
	print("#                                                                  #")
	print("#  Python terminal can be found in "+pythonCommand+"               #")
	print("#  To run BioSANS GUI, type either one of the following commands:  #")
	print("#                                                                  #")
	print("#  BioSANS                                                         #")
	print("#  "+pythonCommand+" -m BioSANS2020                                #")
	print("#  "+pythonCommand+" -m BioSANS2020.BioSANS                        #")
	print("#                                                                  #")
	print("#  To run BioSSL CLI, type either one of the following commands:   #")
	print("#                                                                  #")
	print("#  BioSSL                                                          #")
	print("#  "+pythonCommand+" -m BioSANS2020.BioSSL                         #")
	print("#                                                                  #")
	print("####################################################################")

	Drive = os.path.abspath(os.environ["HOMEPATH"]).split(":")[0]+":"
	try:
		os.mkdir(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS").replace(Drive,"C:"),0o777)
	except:
		pass
	Shcut = open(os.path.abspath(os.environ["HOMEPATH"]+"/BioSANS/BioSANSb.bat").replace(Drive,"C:"),"w")
	Shcut.write("@echo off\n")
	Shcut.write("start "+pythonCommand+"w -m BioSANS2020.BioSANS")
	Shcut.close()
	print("#####################################################")
	print("#                                                   #")
	print("#  You may also launched BioSANS by double click:   #")
	print("#                                                   #")
	print("#  BioSANS shortcut in Desktop                      #")
	print("#  BioSANS shortcut in Start Menu                   #")	
	print("#                                                   #")
	print("#####################################################")
	
	p = Popen("short_cut_create.bat", cwd=os.getcwd())
	stdout, stderr = p.communicate()
		
else:
	print("Please install python 3.7 or higher first")
	
end_na = "None"
while end_na.lower().strip() != "y" and end_na.lower().strip() != "Y":
	try:
		end_na = input("Exit now? Y/N")
	except:
		end_na = raw_input("Exit now? Y/N")