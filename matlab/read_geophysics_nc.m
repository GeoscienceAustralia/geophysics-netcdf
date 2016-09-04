clc;
clear all

ncdir = 'Z:\projects\geophysics_netcdf\ncfiles\';
ncfile = [ncdir 'PIRSA_P702_SAEI_C1_tmi.nc'];

ncfileid = netcdf.open(ncfile,'NC_NOWRITE');
varid = netcdf.inqvarid(ncfileid,'_index_line');
v=netcdf.getvar(ncfileid,varid);

netcdf.close(ncfileid);

