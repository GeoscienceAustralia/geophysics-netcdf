function data = get_line(fileid,vid,start,count,stride)
    
    if(nargin==5)                
        nread = ceil(count/stride);
    else
        stride=1;
        nread=count;
    end  
    
    data = netcdf.getVar(fileid,vid,start,nread,stride);
    fv = get_fill_value(fileid,vid);
    data(data==fv)=NaN;
end
