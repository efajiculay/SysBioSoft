import sys
import os
from subprocess import Popen
import wget

url = "https://www.python.org/ftp/python/3.9.5/python-3.9.5-amd64.exe"
wget.download(url, os.getcwd())
p = Popen("BioSANS_installer.bat", cwd=os.getcwd())
stdout, stderr = p.communicate()