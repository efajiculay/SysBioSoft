import sys
from os import getcwd, path, environ, system
from subprocess import Popen, PIPE, check_output
from wget import download

try:
    gcwd = sys._MEIPASS
except:
    gcwd = getcwd()

def install_BioSANS_old():
    Bcwd = path.join(gcwd,"BioSANSLib")
    file_path = path.join(gcwd,"BioSANSLib","BioSANS_installer.sh")

    print("\n\nInstalling python and necessary libraries\n\n")
    p = Popen(['bash',file_path], cwd=Bcwd)
    stdout, stderr = p.communicate()

def install_BioSANS():
    commands = [
        'xcode-select —install',
        '/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"',
        'brew install python',
        'python3 -m pip install --upgrade pip',
        'pip3 install pillow',
        'pip3 install numpy',
        'pip3 install matplotlib==3.3.3',
        'pip3 install simps',
        'pip3 install spicy',
        'pip3 install sdeint',
        'pip3 install func_timeout',
        'pip3 install python-libsbml',
		'pip3 install pandas',
		'pip3 install applescript',
        'Brew install python-tk',
        'python3 -m pip install -i https://test.pypi.org/simple/ BioSANS2020-efajiculay'
    ]

    for command in commands:
        system(command)
    

def generate_shortcut():
    A = check_output(['pip3', 'show', 'BioSANS2020-efajiculay'])
    A = str(A).split("\\n")

    install_dir = ""
    for row in A:
        line = row.split(":")
        if line[0].strip() == "Location":
            install_dir = ("".join(line[1:]).strip("\\r\\n").replace("c\\","c:/").replace("\\","/").replace("//","/"))
    print("BioSANS installation directory is :",install_dir)

    if install_dir != "":
        install_dir = str(install_dir)
        print()
        print("To run BioSANS GUI, type the following commands")
        print()
        print("python3"+" "+install_dir+"/BioSANS2020/BioSANS.py")
        print()
        print()
        print("To run BioSANS SSL console, type the following commands")
        print()
        print("python3"+" "+install_dir+"/BioSANS2020/BioSSL.py")
        print()
        
        shcutPath = "/Users/"+environ['USER']
        Shcut = open(shcutPath+"/Desktop/BioSANS.sh","w")
        Shcut.write("python3 "+install_dir+"/BioSANS2020/BioSANS.py")
        Shcut.close()
        print()
        print("You may also launched BioSANS by double clicking the BioSANS in the desktop and in applications")
        print()

        system("cp "+ gcwd +"/BioSANSLib/BioSANS "+shcutPath+"/Desktop/")
        system("cp "+ gcwd +"/BioSANSLib/BioSANS /Applications/")
    else:
        print("Installation failed or directory of installation not properly grabbed")


#install_BioSANS_old()
install_BioSANS()
generate_shortcut()
