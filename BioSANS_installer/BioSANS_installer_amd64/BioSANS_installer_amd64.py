import sys
from os import getcwd, path
from subprocess import Popen
from wget import download

Bcwd = path.join(getcwd(),"BioSANSLib")
file_path = path.join(getcwd(),"BioSANSLib","BioSANS_installer.bat")
print(file_path)

print("\nDownloading python from repository\n\n")
url = "https://www.python.org/ftp/python/3.9.5/python-3.9.5-amd64.exe"
download(url, Bcwd)

print("\n\nInstalling python and necessary libraries\n\n")
p = Popen(file_path, cwd=Bcwd)
stdout, stderr = p.communicate()