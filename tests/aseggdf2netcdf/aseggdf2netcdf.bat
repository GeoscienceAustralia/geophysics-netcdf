@echo off

REM set path=C:\Microsoft_HPC_Pack_2012\Bin;%path%
REM set path=..\..\..\..\fftw3.2.2.dlls\64bit;%path%
REM set path=..\..\..\bin\x64\Release\;%path%
set path=Z:\code\repos\geophysics-netcdf\bin\x64\Release;%path%

aseggdf2netcdf.exe aseggdf2netcdf.con

pause

