function data = get_by_stdname(fileid,stdname)
vid = get_vid_by_stdname(fileid,stdname);
data = netcdf.getVar(fileid,vid);
fv = get_fill_value(fileid,vid);
data(data==fv)=NaN;
end
