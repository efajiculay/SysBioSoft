# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(['MacOSX_installer.py'],
             pathex=['/Users/erickson/Desktop/BioSANS_Mac_Installer'],
             binaries=[],
             datas=[('BioSANSLib/BioSANS.icns', 'BioSANSLib'), ('BioSANSLib/BioSANS.zip', 'BioSANSLib'), ('BioSANSLib/BioSANS_installer.sh', 'BioSANSLib'), ('BioSANSLib/BioSANS.py', 'BioSANSLib'), ('BioSANSLib/BioSANS.sh', 'BioSANSLib')],
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
          name='MacOSX_installer',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True , icon='BioSANSLib/BioSANS.icns')
