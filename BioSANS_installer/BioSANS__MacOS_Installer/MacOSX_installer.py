import sys
from os import getcwd, environ, system
import zipfile
import subprocess

try:
    gcwd = sys._MEIPASS
except:
    gcwd = getcwd()



pycom = '$HOME/BioSANS2022/miniconda/bin/python3'

def install_BioSANS():

    commands = [
        'xcode-select --install',
        '/bin/bash -c "$(curl -fsSL https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o ~/miniconda.sh)"',
        'bash ~/miniconda.sh -b -p $HOME/BioSANS2022/miniconda',
        pycom+' -m pip install --upgrade pip',
        pycom+' -m pip install pillow',
        pycom+' -m pip install numpy',
        pycom+' -m pip install matplotlib==3.3.3',
        pycom+' -m pip install simps',
        pycom+' -m pip install spicy',
        pycom+' -m pip install func_timeout',
        pycom+' -m pip install python-libsbml',
        pycom+' -m pip install pandas',
        pycom+' -m pip install applescript',
        pycom+' -m pip install BioSANS2020 --upgrade'
    ]

    for command in commands:
        system(command)


def generate_shortcut():

    print("")
    print("To run BioSANS GUI, type the following commands in terminal")
    print("")
    print(pycom+" -m BioSANS2020")
    print("")
    print("")
    print("To run BioSANS SSL console, type the following commands")
    print("")
    print(pycom+" -m BioSANS2020.BioSSL")
    print("")
    
    shcutPath = "/Users/"+environ['USER']
    Shcut = open(shcutPath+"/Desktop/BioSANS.sh","w")
    Shcut.write(pycom+" -m BioSANS2020")
    Shcut.close()
    print("")
    print("You may also launched BioSANS by double clicking the BioSANS in the desktop and in applications")
    print("")

    with zipfile.ZipFile(gcwd +"/BioSANSLib/BioSANS.zip", 'r') as zip_ref:
        zip_ref.extractall(shcutPath+"/Desktop/")
        subprocess.call(['chmod', '-R', '777', shcutPath+"/Desktop/BioSANS.app"])
        subprocess.call(['rm','-r',shcutPath+"/Desktop/__MACOSX"])

    with zipfile.ZipFile(gcwd +"/BioSANSLib/BioSANS.zip", 'r') as zip_ref:
        zip_ref.extractall("/Applications/")
        subprocess.call(['chmod', '-R', '777', "/Applications/BioSANS.app"])
        subprocess.call(['rm','-r',"/Applications/__MACOSX"])

    # system("cp "+gcwd +"/BioSANSLib/BioSANS "+shcutPath+"/Desktop/")
    # system("cp "+gcwd +"/BioSANSLib/BioSANS /Applications/")


install_BioSANS()
generate_shortcut()
