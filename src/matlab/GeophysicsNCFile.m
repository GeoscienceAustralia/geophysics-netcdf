classdef GeophysicsNCFile < handle
    
    properties (SetAccess = private)
        path;
        ncid;
        chunksize;
        info;
        start;
        count;
    end
    
    methods
        
        function F = GeophysicsNCFile(filepath,mode,chunksize)
            if(nargin<2); mode = 'NOWRITE'; end;
            if(nargin<3); chunksize = 5e6; end;
            [F.chunksize,F.ncid] = netcdf.open(filepath,mode,chunksize);
            F.path = filepath;
            F.info = ncinfo(F.path);
        end
        
        function close(F)
            netcdf.close(F.ncid);
        end
        
        function ncdisp(F)
            ncdisp(F.path);
        end
        
        function I = ncinfo(F)
            I = ncinfo(F.path);
        end
        
        function init_start_count(F)
            line_index = get_by_name(F.ncid,'line_index');           
            [u,F.start] = unique(line_index,'stable');
            F.start = [F.start-1]';    
            F.count = diff([F.start length(line_index)]);            
        end
        
        function check_start_count(F)
            if(isempty(F.start))
                init_start_count(F);
            end            
        end
        
        function n = nlines(F)
            F.check_start_count();
            n = length(F.start);
        end
        
        function n = nsamples(F,lineindex)
            F.check_start_count();
            n = F.count(lineindex);
        end

        function list_dims(F)
            [ndims, nvars, nglobalatts, unlimdimID] = netcdf.inq(F.ncid);
            for di=0:1:ndims-1
                [dimname, dimlen] = netcdf.inqDim(F.ncid,di);
                disp(sprintf('%s: %d',dimname,dimlen));
            end
        end
        
        function list_vars(F)
            [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(F.ncid);
            for vi=0:1:numvars-1
                [varname, xtype, dimids, numatts] = netcdf.inqVar(F.ncid,vi);
                disp(varname);
            end
        end
        
        function s = dim_string_name(F,dimids)
            s='(';
            for di=1:1:length(dimids)
                [dimname, dimlen] = netcdf.inqDim(F.ncid,dimids(di));
                s=[s sprintf('%s',dimname)];
                if(di<length(dimids));
                    s=[s ' x '];
                end
            end
            s=[s ')'];
        end
        
        function s = dim_string_len(F,dimids)
            s='(';
            for di=1:1:length(dimids)
                [dimname, dimlen] = netcdf.inqDim(F.ncid,dimids(di));
                s=[s sprintf('%d',dimlen)];
                if(di<length(dimids));
                    s=[s ' x '];
                end
            end
            s=[s ')'];
        end
        
        function s = dim_string(F,dimids)
            s = [dim_string_name(F,dimids) ':' dim_string_len(F,dimids)];
        end
        
        function list_vars_ex(F)
            [numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(F.ncid);
            for vi=0:1:numvars-1
                [varname, xtype, dimids, numatts] = netcdf.inqVar(F.ncid,vi);
                ty = F.nc_type_string(xtype);
                s = sprintf('%s %s %s',ty,varname,F.dim_string(dimids));
                disp(s);
            end
        end                
        
        function var = get_var_by_name(F,varname)
            var = [];
            for i=1:1:length(F.info.Variables)                
                if(strcmp(varname,F.info.Variables(i).Name))
                    var = F.info.Variables(i);
                    break;
                end
            end                                    
        end        
        
        function var = get_var_by_attval(F,attname,attval)
            var = [];
            for i=1:1:length(F.info.Variables)                
                att = F.get_varatt(F.info.Variables(i),attname);
                if(~isempty(att))
                    if(strcmp(att.Value,attval))
                        var = F.info.Variables(i);
                        break;
                    end
                end
            end                                    
        end        
        
        function var = get_var_by_stdname(F,stdname)            
            var = get_var_by_attval(F,'standard_name',stdname);
        end        
        
        function var = get_var_by_longname(F,longname)            
            var = get_var_by_attval(F,'long_name',longname);
        end  
        
        function d = get_vardata(F,var)
            varid = netcdf.inqVarID(F.ncid,var.Name);
            if(varid<0)
                d = [];  
            else
                d = netcdf.getVar(F.ncid,varid);
                fv = F.get_fill_value(var);
                if(~isempty(fv))
                    d(d==fv)=NaN;
                end                            
            end
        end
        
        function d = get_vardata_by_name(F,varname)
            var = F.get_var_by_name(varname);
            d = F.get_vardata(var);
        end
        
        function d = get_vardata_lineindex(F,var,lineindex,stride)
            
            varid = netcdf.inqVarID(F.ncid,var.Name);
            if(varid < 0)
                msg = sprintf('Invalid varid for %s',var.Name);
                error(msg);
            end
                
            if(nargin==4)                
                numtoread = ceil(F.count(lineindex)/stride);
            else
                stride=1;
                numtoread = F.count(lineindex);
            end  

            d = netcdf.getVar(F.ncid,varid,F.start(lineindex),numtoread,stride);            
            fv = F.get_fill_value(var);
            if(~isempty(fv))
                d(d==fv)=NaN;
            end
        end

        function bp = get_bounding_polygon(F)
            var = F.get_var_by_name('bounding_polygon');
            bp  = F.get_vardata(var);                        
        end
        
        function [xstart,xend,ystart,yend] = get_line_start_end_points(F)            
            xstart = F.get_vardata_by_name('longitude_first');
            xend   = F.get_vardata_by_name('longitude_last');
            ystart = F.get_vardata_by_name('latitude_first');
            yend   = F.get_vardata_by_name('latitude_last');
        end
        
        
        
    end
    methods (Static)
        function s = nc_type_string(xtype)
            ty = {'int8', 'char','short','int32', 'float','double','uint8','ushort','uint'};
            s = char(ty(xtype));
        end
        
        function att = get_varatt(var,attname)
            att = [];
            for i=1:1:length(var.Attributes)                            
                  if(strcmp(var.Attributes(i).Name,attname))
                     att = var.Attributes(i);
                     break;
                  end
            end                        
        end
        
        function fv = get_fill_value(var)
            fv=[];
            att = GeophysicsNCFile.get_varatt(var,'_FillValue');            
            if(~isempty(att));
                fv = att.Value;
            end
        end
        
    end
end




