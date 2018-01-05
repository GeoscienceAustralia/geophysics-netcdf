function attval = get_att_value(fileid,vid,attname)
[varname, xtype, dimids, numatts] = netcdf.inqVar(fileid,vid);
for ai=0:1:numatts-1
    name = netcdf.inqAttName(fileid,vid,ai);
    if(strcmp(name,attname))
        attval = netcdf.getAtt(fileid,vid,attname);
        return;
    end
end
attval = [];
end