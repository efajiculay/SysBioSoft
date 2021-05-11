import sys
import os
import subprocess
import pkg_resources

try:
	os.system("sudo apt update --user")
except:
	pass
try:
	os.system("sudo apt install python3-pip")
except:
	s.system("sudo apt install python-pip")
finally:
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
		os.system("pip3 install "+item+" --user")
	except:
		os.system("pip install "+item+" --user")
	finally:
		pass

try:
	os.system("sudo apt-get install python3-tk"+" --user")
except:
	pass
try:
	os.system("sudo apt-get install python-tk"+" --user")
except:
	pass
	
try:
	os.system("sudo apt-get install python3-pil python3-pil.imaget")
except:
	pass
	
try:
	os.system("pip3 install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay")
except:
	os.system("pip install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay")
finally:
	pass	

try:
	A = subprocess.check_output(['pip3', 'show', 'BioSANS2020-efajiculay'])
except subprocess.CalledProcessError as err:
	A = subprocess.check_output(['pip', 'show', 'BioSANS2020-efajiculay'])
finally:
	pass
	
    
A = str(A).split("\\n")
install_dir = ""
for row in A:
	line = row.split(":")
	if line[0].strip() == "Location":
		install_dir = line[1].strip()

		
if install_dir != "":
	install_dir = str(install_dir)
	print()
	print("To run BioSANS GUI, type the following commands")
	print()
	print("python3 "+install_dir+"/BioSANS2020/BioSANS.py")
	print()
	print("or")
	print()
	print("python "+install_dir+"/BioSANS2020/BioSANS.py")
	print()
	print("To run BioSANS SSL console, type the following commands")
	print()	
	print("python3 "+install_dir+"/BioSANS2020/BioSSL.py")
	print()
	print("or")
	print()
	print("python "+install_dir+"/BioSANS2020/BioSSL.py")
	print()
else:
	print("software installed but directory of installation not automatically grabbed")
	print("Locate installation directory using\n")
										 
			
	print("pip3 show BioSANS2020-efajiculay")
	print("or")
	print("pip show BioSANS2020-efajiculay")
	print()
	print("To run BioSANS GUI, type the following commands")
	print()
	print("python3 "+install_dir+"/BioSANS2020/BioSANS.py")
	print()
	print("or")
	print()
	print("python "+install_dir+"/BioSANS2020/BioSANS.py")
	print()
	print("To run BioSANS SSL console, type the following commands")
	print()	
	print("python3 "+install_dir+"/BioSANS2020/BioSSL.py")
	print()
	print("or")
	print()
	print("python "+install_dir+"/BioSANS2020/BioSSL.py")
	print()