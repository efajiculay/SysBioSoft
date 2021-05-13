:: https://silentinstallhq.com/python-3-9-silent-install-how-to-guide/
:: https://www.tutorialspoint.com/batch_script/
:: @echo off

start python-3.9.5.exe /uninstall /quiet
mkdir C:\BioSANS2021b
start /WAIT .\python-3.9.5.exe /quiet .\unattend.xml
start /B /WAIT C:\BioSANS2021b\Python39\python.exe Window_Installer.py
