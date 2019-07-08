function vid = get_vid_by_attval(fileid,attname,attval)
    [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(fileid);
    vid = -1;
    for vi=0:1:numvars-1
        [varname, xtype, dimids, numatts] = netcdf.inqVar(fileid,vi);
        val = get_att_value(fileid,vi,attname);
        if(strcmp(val,attval))
            vid = vi;
            return;
        end
    end
end

