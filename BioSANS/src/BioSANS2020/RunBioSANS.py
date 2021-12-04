from os import system
from subprocess import Popen
from sys import executable


def BioSANS():
    Popen([str(executable), "-m", "BioSANS2020.BioSANS"])


def BioSSL():
    system(str(executable)+" -m BioSANS2020.BioSSL")
