
"""
import sys
from subprocess import check_output
from os import system

A = check_output(['pip3', 'show', 'BioSANS2020-efajiculay'])
#print(A)
A = str(A).split("\\n")
#print(A)

install_dir = ""
for row in A:
    line = row.split(":")
    if line[0].strip() == "Location":
        install_dir = ("".join(line[1:]).strip("\\r\\n").replace("c\\","c:/").replace("\\","/").replace("//","/"))

if install_dir != "":
    install_dir = str(install_dir)
    system("python3 "+install_dir+"/BioSANS2020/BioSANS.py")
"""

from os import system


pycom = '$HOME/BioSANS2022/miniconda/bin/python3'
system(pycom+" -m BioSANS2020")

