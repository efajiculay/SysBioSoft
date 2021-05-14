:: start /B pyinstaller ^
:: --onefile ^
:: --icon BioSANS.ico ^
:: --add-data BioSANS.ico;. ^
:: --add-data BioSANS_installer.bat;. ^
:: --add-data short_cut_create.bat;. ^
:: --add-data unattend.xml;. ^
:: --add-data Window_Installer.py;. ^
:: --clean BioSANS_installer_amd64.py

start /B pyinstaller --onefile --icon BioSANSLib/BioSANS.ico --clean BioSANS_installer_amd64.py