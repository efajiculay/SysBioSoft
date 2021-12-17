"""

                  This module is the RunBioSANS module

This can be used to run BioSANS and BioSSL byb running the folloeing.

BioSSL()
BioSANS()


"""

from os import system
from subprocess import Popen
from sys import executable


def BioSANS():
    """This function launched BioSANS"""
    Popen([str(executable), "-m", "BioSANS2020.BioSANS"])


def BioSSL():
    """This function launched BioSSL"""
    system(str(executable) + " -m BioSANS2020.BioSSL")
