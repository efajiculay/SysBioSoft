@echo off

xcopy %cd%\BioSANS.ico %USERPROFILE%\BioSANS

set SCRIPT="%TEMP%\%RANDOM%-%RANDOM%-%RANDOM%-%RANDOM%.vbs"

echo Set oWS = WScript.CreateObject("WScript.Shell") >> %SCRIPT%

echo sLinkFile = "%USERPROFILE%\Desktop\BioSANSb.lnk" >> %SCRIPT%
echo Set oLink = oWS.CreateShortcut(sLinkFile) >> %SCRIPT%
echo oLink.IconLocation = "%USERPROFILE%\BioSANS\BioSANS.ico" >> %SCRIPT%
echo oLink.TargetPath = "%USERPROFILE%\BioSANS\BioSANSb.bat" >> %SCRIPT%
echo oLink.WorkingDirectory = "%USERPROFILE%\BioSANS" >> %SCRIPT% 
echo oLink.Save >> %SCRIPT%

echo sLinkFile = "%AppData%\Microsoft\Windows\Start Menu\Programs\BioSANSb.lnk" >> %SCRIPT%
echo Set oLink = oWS.CreateShortcut(sLinkFile) >> %SCRIPT%
echo oLink.IconLocation = "%USERPROFILE%\BioSANS\BioSANS.ico" >> %SCRIPT%
echo oLink.TargetPath = "%USERPROFILE%\BioSANS\BioSANSb.bat" >> %SCRIPT%
echo oLink.WorkingDirectory = "%USERPROFILE%\BioSANS" >> %SCRIPT% 
echo oLink.Save >> %SCRIPT%

cscript /nologo %SCRIPT%
del %SCRIPT%


