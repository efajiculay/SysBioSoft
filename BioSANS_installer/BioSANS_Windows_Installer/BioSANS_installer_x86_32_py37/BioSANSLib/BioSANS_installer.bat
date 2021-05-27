:: https://silentinstallhq.com/python-3-9-silent-install-how-to-guide/
:: https://www.tutorialspoint.com/batch_script/
:: @echo off

:: start python-3.7.4.exe /uninstall /quiet
mkdir C:\BioSANS2021\Python37x86_32
start /WAIT .\python-3.7.4.exe /quiet .\unattend.xml
start /B /WAIT C:\BioSANS2021\Python37x86_32\Python37\python.exe Window_Installer.py
