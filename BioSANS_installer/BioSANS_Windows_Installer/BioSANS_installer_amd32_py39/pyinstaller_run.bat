start /B pyinstaller ^
--onefile ^
--icon BioSANSLib/BioSANS.ico ^
--add-data BioSANSLib/BioSANS.ico;BioSANSLib ^
--add-data BioSANSLib/BioSANS_installer.bat;BioSANSLib ^
--add-data BioSANSLib/short_cut_create.bat;BioSANSLib ^
--add-data BioSANSLib/unattend.xml;BioSANSLib ^
--add-data BioSANSLib/Window_Installer.py;BioSANSLib ^
--clean BioSANS_installer_amd32.py

:: start /B pyinstaller --onefile --icon BioSANSLib/BioSANS.ico --clean BioSANS_installer_amd32.py