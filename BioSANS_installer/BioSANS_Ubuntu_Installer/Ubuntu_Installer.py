import sys
from os import getcwd, path, environ, system
import subprocess
import pkg_resources

try:
    gcwd = sys._MEIPASS
except:
    gcwd = getcwd()


try:
	system("sudo apt update")
except:
	pass
try:
	system("sudo apt install python3-pip")
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
		system("pip3 install "+item+" --user")
	except:
		system("pip install "+item+" --user")
	finally:
		pass

try:
	system("sudo apt-get install python3-tk"+" --user")
except:
	pass
try:
	system("sudo apt-get install python-tk"+" --user")
except:
	pass
	
try:
	system("sudo apt-get install python3-pil python3-pil.imaget")
except:
	pass
	
try:
	system("pip3 install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay")
except:
	system("pip install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay")
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
	print("To run BioSANS SSL console, type the following commands")
	print()	
	print("python3 "+install_dir+"/BioSANS2020/BioSSL.py")
	print()
	
	shcutPath = environ['HOME']
	Shcut = open(shcutPath+"/Desktop/BioSANS.sh","w")
	Shcut.write("python3 "+install_dir+"/BioSANS2020/BioSANS.py")
	Shcut.close()
	print()
	print("You may also launched BioSANS by double clicking the BioSANS shortcut in the desktop and in applications")
	print()

	system("mkdir "+shcutPath+"/BioSANSLib")
	system("cp "+ gcwd +"/BioSANSLib/BioSANS "+shcutPath+"/BioSANSLib")
	system("cp "+ gcwd +"/BioSANSLib/BioSANS.xpm "+shcutPath+"/BioSANSLib")
		
	Shcut = open(shcutPath+"/Desktop/BioSANS.desktop","w")
	Shcut.write("[Desktop Entry]\n")
	Shcut.write("Name=BioSANS\n")
	Shcut.write("Comment=BioSANS2020\n")
	Shcut.write("Exec="+shcutPath+"/BioSANSLib/BioSANS"+"\n")
	Shcut.write("Icon="+shcutPath+"/BioSANSLib/BioSANS.xpm"+"\n")
	Shcut.write("Terminal=false\n")
	Shcut.write("Type=Application\n")
	Shcut.write("Categories=Development;\n")
	Shcut.write("StartupNotify=true\n")
	Shcut.write("NoDisplay=true\n")
	Shcut.close()	
		
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
