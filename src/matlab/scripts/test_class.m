clc;
clear all;

addpath('Z:\code\repos\geophysics-netcdf\src\matlab\');
ncfilepath = 'Z:\projects\geophysics_netcdf\ncfiles\v2\P441MAG.nc';

F = GeophysicsNCFile(ncfilepath);
%F.disp();
%i = F.info();
%F.list_dims();
%F.list_vars_ex();
var = F.get_var_by_name('fiducial');
vid = F.get_varid_by_name('fiducial');
F.close();

