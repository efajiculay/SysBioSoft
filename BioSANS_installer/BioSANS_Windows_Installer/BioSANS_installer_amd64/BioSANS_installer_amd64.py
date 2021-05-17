import sys
from os import getcwd, path
from subprocess import Popen, PIPE
from wget import download

gcwd = sys._MEIPASS 
#gcwd = getcwd()

Bcwd = path.join(gcwd,"BioSANSLib")
file_path = path.join(gcwd,"BioSANSLib","BioSANS_installer.bat")

url = "https://www.python.org/ftp/python/3.9.5/python-3.9.5-amd64.exe"
if not path.isfile(path.join(Bcwd,"python-3.9.5-amd64.exe")):
	print("\nDownloading python from repository\n\n")
	download(url, Bcwd)

print("\n\nInstalling python and necessary libraries\n\n")
p = Popen(file_path, cwd=Bcwd)
stdout, stderr = p.communicate()