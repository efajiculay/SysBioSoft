# -*- mode: python ; coding: utf-8 -*-


block_cipher = None


a = Analysis(['BioSANS_win10_amd32_py39.py'],
             pathex=['C:\\Users\\Erickson\\Documents\\GitHub\\SysBioSoft\\BioSANS_installer\\BioSANS_Windows_New_Installer\\BioSANS_win10_amd32_py39'],
             binaries=[],
             datas=[],
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
          name='BioSANS_win10_amd32_py39',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          upx_exclude=[],
          runtime_tmpdir=None,
          console=True , icon='BioSANS.ico')
