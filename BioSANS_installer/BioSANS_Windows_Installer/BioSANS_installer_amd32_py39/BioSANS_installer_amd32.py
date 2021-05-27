import sys
from subprocess import Popen, PIPE, STDOUT, REALTIME_PRIORITY_CLASS
from os import getcwd, path, system, mkdir
from wget import download

try:
	gcwd = sys._MEIPASS
except:	
	gcwd = getcwd()

Bcwd     = path.join(gcwd,"BioSANSLib")
pyInsexe = path.join(Bcwd,"python-3.9.5.exe")
xmlset   = path.join(Bcwd,"unattend.xml")
WinIns   = path.join(Bcwd,"Window_Installer.py")
shcut    = path.join(Bcwd,"short_cut_create.bat")

url = "https://www.python.org/ftp/python/3.9.5/python-3.9.5.exe"
if not path.isfile(path.join(Bcwd,"python-3.9.5.exe")):
	print("\nDownloading python from repository\n\n")
	download(url, Bcwd)

print("\n\nInstalling python and necessary libraries\n\n")

try:
	mkdir("C:/BioSANS2021/python39amd32",0o777)
except:
	pass
#system(pyInsexe+" /quiet "+xmlset)

p = Popen([pyInsexe, "/quiet", xmlset], stdout = PIPE, stderr = STDOUT, shell = True, creationflags=REALTIME_PRIORITY_CLASS)

done = False
count = 1
while not done:
	try:
		outs, errs = p.communicate(timeout=2)
		done = True
	except:
		print("processing isolated local python installation : Approximately",min(100,round(100*count/45,0)),"%")
	count = count + 1
print("processing isolated local python installation : Approximately",100,"%")

print("\nBioSANS Installation started :\n")
#system("C:/BioSANS2021/python39amd32/Python39/python.exe "+WinIns)
p = Popen(["C:/BioSANS2021/python39amd32/Python39/python.exe", WinIns], stdout = PIPE, stderr = STDOUT, shell = True, creationflags = REALTIME_PRIORITY_CLASS)

done = False
count = 1
while not done:
	try:
		outs, errs = p.communicate(timeout=2)
		done = True
	except:
		print("processing BioSANS installation : Approximately",min(100,round(100*count/70,0)),"%")
	count = count + 1
print("processing BioSANS installation : Approximately",100,"%\n")

print()
print("##################################################################################")
print("#                                                                                #")
print("#  Python terminal can be found in C:/BioSANS2021/python39amd32/Python39/python  #")
print("#  To run BioSANS GUI, type either one of the following 3 commands:              #")
print("#                                                                                #")
print("#  BioSANS                                                                       #")
print("#  C:/BioSANS2021/python39amd32/Python39/python -m BioSANS2020                   #")
print("#  C:/BioSANS2021/python39amd32/Python39/python -m BioSANS2020.BioSANS           #")
print("#                                                                                #")
print("#  To run BioSSL CLI, type either one of the following 2 commands:               #")
print("#                                                                                #")
print("#  BioSSL                                                                        #")
print("#  C:/BioSANS2021/python39amd32/Python39/python -m BioSANS2020.BioSSL            #")
print("#                                                                                #")
print("##################################################################################")
print()

p = Popen(shcut, cwd=Bcwd)
stdout, stderr = p.communicate()

print()
print("#####################################################")
print("#                                                   #")
print("#  You may also launched BioSANS by double click:   #")
print("#                                                   #")
print("#  BioSANS shortcut in Desktop                      #")
print("#  BioSANS shortcut in Start Menu                   #")	
print("#                                                   #")
print("#####################################################")
print()

end_na = "None"
while end_na.lower().strip() != "y" and end_na.lower().strip() != "Y":
	try:
		end_na = input("Exit now? Y/N")
	except:
		end_na = raw_input("Exit now? Y/N")
