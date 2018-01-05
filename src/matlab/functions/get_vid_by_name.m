function vid = get_vid_by_name(fileid,vname)
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(fileid);
for vi=0:1:numvars-1
    [varname, xtype, dimids, numatts] = netcdf.inqVar(fileid,vi);
    if(strcmp(varname,vname))
        vid = vi;
        return;
    end
end
vid = -1;
end
