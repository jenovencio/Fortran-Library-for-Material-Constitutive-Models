@echo off
rem generateUserMatLib.bat (list of subroutines)
rem generateUserMatLib.bat libfiles.txt
if exist usermat.f (del usermat.f)

for /f %%G in (%1) do ( type %%G >> usermat.f & echo c> fake.txt & type fake.txt >> usermat.f & del fake.txt)
rem findstr /v "!end code" usermat.f > usermatlib2.f90
rem del usermat.f
rem move usermatlib2.f90 usermat.f

echo ---------------------------
echo usermat.f was created
echo ---------------------------
