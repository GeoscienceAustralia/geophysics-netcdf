function fv = get_fill_value(fileid,vid)
[varname, xtype, dimids, numatts] = netcdf.inqVar(fileid,vid);
for ai=0:1:numatts-1
    attname = netcdf.inqAttName(fileid,vid,ai);
    if(strcmp(attname,'_FillValue'))
        fv = netcdf.getAtt(fileid,vid,attname);
        return
    end
end
fv=NaN;
end
