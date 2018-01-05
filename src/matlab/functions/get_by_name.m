function data = get_by_name(fileid,vname)
vid = get_vid_by_name(fileid,vname);
data = netcdf.getVar(fileid,vid);
fv = get_fill_value(fileid,vid);
data(data==fv)=NaN;
end
