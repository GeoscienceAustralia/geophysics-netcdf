function data = get_by_name(fileid,vname)
    vid = get_vid_by_name(fileid,vname);
    if(vid<0)
        data = [];        
    else
        data = netcdf.getVar(fileid,vid);
        fv = get_fill_value(fileid,vid);
        data(data==fv)=NaN;
    end
end
