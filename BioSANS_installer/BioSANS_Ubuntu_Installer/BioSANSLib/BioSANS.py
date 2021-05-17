import sys
from subprocess import check_output
from os import system


#gio set Your_desktop_file.desktop "metadata::trusted" yes

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
