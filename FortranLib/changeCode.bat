@echo off
rem changeCode.bat oldstr newstr list
rem set oldstr = real
rem newstr = DOUBLE PRECISION
rem list = list.txt
rem generateUserMatLib.bat %oldstr% %newstr% %list%
if exist usermat.f (del usermat.f)

for /f %%G in (%3) do ( BatchSubstitute.bat %1 %2 %%G)

echo ---------------------------
echo Code was sucessfull created
echo ---------------------------