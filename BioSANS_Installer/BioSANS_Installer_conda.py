import sys
import os
import subprocess
import pkg_resources

#print(os.environ["HOMEPATH"])
#print(sys.version)
#conda env remove -n BioSANS
#conda create -n BioSANS python==3.8.8
#conda activate BioSANS

python_ok = False
if sys.version_info[0] == 3 and sys.version_info[1] >= 7:
	python_ok = True

pythonCommand = "python"
pipCommand = "pip"

if python_ok:

	try:
		os.system(pipCommand+" "+"install --upgrade pip --user")
	except:
		pass

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
			os.system("conda install "+item+"")
		finally:
			os.system(pipCommand+" "+"install "+item+"")		

	try:
		os.system(pipCommand+" "+"install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay --user")
	finally:
		pass

	A = subprocess.check_output([pipCommand, 'show', 'BioSANS2020-efajiculay'])
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

		Shcut = open(os.environ["HOMEPATH"]+"/Desktop/BioSANS.bat","w")
		Shcut.write("@echo off\n")
		Shcut.write("start "+sys.executable.replace(".exe","w.exe")+" "+install_dir+"/BioSANS2020/BioSANS.py")
		Shcut.close()
		print()
		print("You may also launched BioSANS by double clicking the BioSANS.bat file in the desktop")
		print()		
	else:
		print("software installed but directory of installation not automatically grabbed")
		print("Locate installation directory using\n")
		print(pipCommand+" "+"show BioSANS2020-efajiculay")
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
		
else:
	print("Please install python 3.7 or higher first")
	
end_na = "None"
while end_na.lower().strip() != "y" and end_na.lower().strip() != "Y":
	try:
		end_na = input("Exit now? Y/N")
	except:
		end_na = raw_input("Exit now? Y/N")