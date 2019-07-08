function [istart icount] = get_start_count(fileid)    
    line_index = get_by_name(fileid,'line_index');           
    [u,istart] = unique(line_index,'stable');
    istart = [istart-1]';    
    icount = diff([istart length(line_index)]);            
end
