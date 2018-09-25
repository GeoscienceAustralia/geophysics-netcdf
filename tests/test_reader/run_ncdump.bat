@echo off


set path=Y:\ops\gap\software\netcdf\netCDF-4.3.3.1\bin;%path%

FOR /f %%F IN ('dir /b *.nc') DO (
  echo %%F
  ncdump.exe -h %%F > %%F.txt  
)


pause