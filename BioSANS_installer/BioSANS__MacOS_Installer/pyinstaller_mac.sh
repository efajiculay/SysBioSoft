pyinstaller \
--onefile \
--icon BioSANSLib/BioSANS.icns \
--add-data BioSANSLib/BioSANS.icns:BioSANSLib \
--add-data BioSANSLib/BioSANS.zip:BioSANSLib \
--add-data BioSANSLib/BioSANS_installer.sh:BioSANSLib \
--add-data BioSANSLib/BioSANS.py:BioSANSLib \
--add-data BioSANSLib/BioSANS.sh:BioSANSLib \
--clean \
MacOSX_installer.py
