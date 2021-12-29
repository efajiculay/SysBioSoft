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
	pass
	
try:
	system("sudo apt install python-pip")
except:
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
	system("pip3 install BioSANS2020 --upgrade")
except:
	pass
	
try:
	system("pip install BioSANS2020 --upgrade")
except:
	pass	

print()
print("#######################################################################")
print("#                                                                     #")
print("#  To run BioSANS GUI, type any of the following 3 commands           #")
print("#                                                                     #")
print("#  python3 -m BioSANS2020.BioSANS                                     #")
print("#  python3 -m BioSANS2020                                             #")
print("#  BioSANS                                                            #")
print("#                                                                     #")
print("#  To run BioSANS SSL console, type any of the following 2 commands   #")
print("#                                                                     #")	
print("#  python3 -m BioSANS2020.BioSSL                                      #")
print("#  BioSSL                                                             #")
print("#                                                                     #")
print("#######################################################################")

shcutPath = environ['HOME']
Shcut = open(shcutPath+"/Desktop/BioSANS.sh","w")
Shcut.write("python3 -m BioSANS2020.BioSANS")
Shcut.close()

print()
print("#################################################")
print("#                                               #")
print("#  You may also launched BioSANS from           #")
print("#  BioSANS shortcut in the desktop              #")
print("#  activate it first by right click and choose  #")
print("#  allow launching                              #")
print("#                                               #")
print("#################################################")

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
