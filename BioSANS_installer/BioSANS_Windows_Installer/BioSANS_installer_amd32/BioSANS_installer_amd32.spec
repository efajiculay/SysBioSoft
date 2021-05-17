# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['BioSANS_installer_amd32.py'],
             pathex=['C:\\Users\\efaji\\Documents\\GitHub\\SysBioSoft\\BioSANS_installer\\BioSANS_installer_amd32'],
             binaries=[],
             datas=[('BioSANSLib/BioSANS.ico', 'BioSANSLib'), ('BioSANSLib/BioSANS_installer.bat', 'BioSANSLib'), ('BioSANSLib/short_cut_create.bat', 'BioSANSLib'), ('BioSANSLib/unattend.xml', 'BioSANSLib'), ('BioSANSLib/Window_Installer.py', 'BioSANSLib')],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          [],
          name='BioSANS_installer_amd32',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True , icon='BioSANSLib\\BioSANS.ico')
