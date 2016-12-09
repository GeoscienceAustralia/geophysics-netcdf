function data = get_line(fileid,vid,start,count)
    data = netcdf.getVar(fileid,vid,start,count);
    fv = get_fill_value(fileid,vid);
    data(data==fv)=NaN;
end
