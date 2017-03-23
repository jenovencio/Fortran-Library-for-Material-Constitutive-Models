@echo off
setlocal ENABLEDELAYEDEXPANSION
cls
set revn=14.0
set INCLUDE=..\..\..\customize\include;..\..\..\..\customize\include;%INCLUDE%

if "%PROCESSOR_ARCHITECTURE%"=="AMD64" (

   set PLATFORM_DIR=winx64
   set IFORT_PLATFORM=intel64
   set BIT_TARGET=64
   set LIB=..\..\lib\winx64;..\..\..\lib\winx64;!LIB!

) else if "%PROCESSOR_ARCHITEW6432%"=="AMD64" (

   set PLATFORM_DIR=winx64
   set IFORT_PLATFORM=intel64
   set BIT_TARGET=64
   set LIB=..\..\lib\winx64;..\..\..\lib\winx64;!LIB!

) else if "%PROCESSOR_ARCHITECTURE%"=="x86" (

   set PLATFORM_DIR=Intel
   set IFORT_PLATFORM=ia32
   set BIT_TARGET=32
   set LIB=..\..\lib\intel;..\..\..\lib\intel;!LIB!
)

set INCLUDE=%AWP_ROOT140%\ansys\customize\Include;%INCLUDE%
set LIB=%AWP_ROOT140%\ansys\Custom\Lib\%PLATFORM_DIR%;%LIB%

call "%IFORT_COMPILER11%\bin\ifortvars.bat" %IFORT_PLATFORM% vs2008

echo.
echo       **************************************************************
echo.
echo       This is the ANSYS Revision %revn% ANSCUST batch file.  It is
echo       used to  link User Programmable Features into  versions of the
echo       ANSYS program on Microsoft Windows %BIT_TARGET%-bit systems.
echo.
echo                   ******   IMPORTANT !!!! ******
echo.
echo       The user subroutine source files provided by ANSYS now reside
echo       in the folder: {InstallDir}\customize\user 
echo.
echo       Please copy the source you wish to modify from this directory
echo       into your working directory [default {InstallDir}\custom\user\%PLATFORM_DIR%]
echo       to include them in your link.  If you wish to use a working directory other
echo       than the default, copy the ANSCUST.BAT, ansysex.def and ansys.lrf files from
echo       {InstallDir}\custom\user\%PLATFORM_DIR%] to your working directory and run
echo       ANSCUST.BAT from there.
echo.
echo       ******************************************************************
echo.
rem pause


:ASKAERO
SET AERO=FALSE
echo.
set ANSW=N
if /I "%ANSW%"=="" GOTO :PLEASE
if /I %ANSW%==Y (
   SET AERO=TRUE
   goto :COMPILE
)

if /I %ANSW%==N (
    goto :COMPILE
)

:PLEASE
echo.
echo "Please answer Y or N!"
echo.
goto :ASKAERO

:COMPILE
if exist *.log del /Q *.log
if exist *.obj del /Q *.obj

if "%BIT_TARGET%"=="64" (

@echo on
for %%p in ( "*.F" ) do ifort -D__EFL -DNOSTDCALL -DARGTRAIL -O3 -DPCWIN64_SYS -DPCWINX64_SYS -DPCWINNT_SYS -DCADOE_ANSYS -Qfpp -nbs -4Yportlib -4Ya /c /Fo.\ /DFORTRAN /MD /W0 %%p >>f_comp.log  2>&1
IF %ERRORLEVEL% GEQ 2 goto COMPILEFAIL

for %%p in ( "*.c*" ) do cl -DNOSTDCALL -DARGTRAIL -Ox -GR -EHsc /DPCWINNT_SYS /DPCWIN64_SYS /DPCWINX64_SYS /DCADOE_ANSYS /DCURVEFIT_EXPORTS /Zi  /c /W3 /D_X86=1 /DOS_WIN32 /DWIN32 /D__STDC__ /MD %%p >>c_comp.log 2>&1
IF %ERRORLEVEL% GEQ 2 goto COMPILEFAIL
@echo off

)

if "%BIT_TARGET%"=="32" (

@echo on
for %%p in ( "*.F" ) do ifort -D__IFL -DNOSTDCALL -DARGTRAIL -O3 -DPCWINNT_SYS -DCADOE_ANSYS -Qfpp -nbs -4Yportlib -4Ya /align:rec4byte /c /Fo.\ /DFORTRAN /MD /W0 %%p >>f_comp.log  2>&1
IF %ERRORLEVEL% GEQ 2 goto COMPILEFAIL

for %%p in ( "*.c*" ) do cl -DNOSTDCALL -DARGTRAIL -Ox -GR -EHsc /DPCWINNT_SYS /DCADOE_ANSYS /DCURVEFIT_EXPORTS /Zi /Zp4 /c /W3 /D_X86=1 /DOS_WIN32 /DWIN32 /D__STDC__ /MD %%p >>c_comp.log 2>&1
IF %ERRORLEVEL% GEQ 2 goto COMPILEFAIL
@echo off

)

:LINK
if exist ANSYS.exe del /Q ANSYS.exe
if "%AERO%" == "TRUE" (

   type ansys.lrf >ansys.lrf.sav
   echo "%AWP_ROOT140%\ansys\Custom\User\%PLATFORM_DIR%\Aeroelastic\*.obj">>ansys.lrf
   @echo on
   link @ansys.lrf
   type ansys.lrf.sav >ansys.lrf
   del /q ansys.lrf.sav

) ELSE (

   @echo on
   link @ansys.lrf

)

@echo off
IF %ERRORLEVEL% GEQ 1 goto LINKFAIL

:: Embed the manifest
@echo on
mt /outputresource:"ANSYS.exe;#1" /manifest "ANSYS.exe.intermediate.manifest"
@echo off
IF %ERRORLEVEL% GEQ 1 goto MANIFESTFAIL

echo.
echo       ************************************************************************
echo.
echo           LINK WAS SUCCESSFUL!
echo.
echo       ************************************************************************
echo.

:COPYDLL
echo.
echo       ****************************************************************************
echo.
echo       THE NEXT QUESTION WILL GIVE YOU THE OPPORTUNITY TO COPY THE NECESSARY
echo       RUNTIME DLLs.  NOTE, THIS ONLY NEEDS TO BE DONE ONCE.
echo       YOU CAN ANSWER "N" FOR ALL SUBSEQUENT INVOCATIONS OF %0.BAT.
echo.
echo       ****************************************************************************
echo.
set ANSW=Y
if /I %ANSW%==Y (
REM   copy /y ..\..\..\bin\%PLATFORM_DIR%\*.dll .
   copy /y "%AWP_ROOT140%\ansys\Bin\%PLATFORM_DIR%\"*.dll .
   goto :END
)
if /I %ANSW%==N (
    goto :END
)
echo.
echo "Please answer Y or N!"
echo.
goto COPYDLL

goto :END

:COMPILEFAIL
@echo off
echo.
echo.       ************************************************************************
echo.
echo.          THERE HAS BEEN A COMPILE FAILURE!  CHECK f_comp.log or c_comp.log
echo.          FOR DETAILS.
echo.
echo.
echo.      ************************************************************************
echo.
goto :END

:LINKFAIL
echo.
echo       ************************************************************************
echo.
echo           THERE HAS BEEN A LINK FAILURE!  IF YOU ARE LINKING DISTRIBUTED ANSYS
echo           VERIFY THE HP-MPI, MS-MPI or INTEL-MPI LIB PATHS ARE THE SAME AS
echo           SPECIFIED IN THE ANSYS.LRF FILE.
echo.
echo.
echo       ************************************************************************
echo.
goto :END

:MANIFESTFAIL
echo.
echo       ****************************************************************************
echo.
echo           THERE HAS BEEN AN ERROR EMBEDDING THE MANIFEST.  NORMALLY THIS IS CAUSED
echo           WHEN THE MICROSOFT MANIFEST TOOL (MT.EXE) CANNOT BE FOUND.  PLEASE CHECK
echo           YOUR PATH IS SET PROPERLY.
echo.
echo       ****************************************************************************
echo.
goto :END


:END

endlocal
