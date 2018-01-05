clc;
clear all;

addpath('Z:\code\repos\geophysics-netcdf\src\matlab');

%ncfile = '..\aem\ncfiles\AUS_10008_WestK_LCI.nc';
%ncfile  = '..\aem\ncfiles\AUS_10009_Ord_Bonaparte_LCI.nc';
ncfile = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AEM_examples/AUS_10008_WestK_LCI.nc';
%ncfile = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AEM_examples/AUS_10009_Ord_Bonaparte_LCI.nc';

fileid = netcdf.open(ncfile,'NOWRITE');
istart = get_by_name(fileid,'index_line');
icount = get_by_name(fileid,'index_count');
nlines = length(istart);
x1  = get_by_name(fileid,'easting_first');
x2  = get_by_name(fileid,'easting_last');
y1  = get_by_name(fileid,'northing_first');
y2  = get_by_name(fileid,'northing_last');

bp     = get_by_name(fileid,'bounding_polygon');
line   = get_by_name(fileid,'line');
nlines = length(istart);

if(true)
    vx = get_vid_by_name(fileid,'easting');
    vy = get_vid_by_name(fileid,'northing');
    
    figure(1);
    maximize_figure();
    tic
    plot([x1 x2]',[y1 y2]',':b');
    hold on;           
    for k=1:1:nlines
        x = get_line(fileid,vx,istart(k),icount(k));
        y = get_line(fileid,vy,istart(k),icount(k));
        plot(x,y,'-g');
    end
    plot(bp(1,:),bp(2,:),'-r.');    
    daspect([1 1 1]);   
    drawnow;
end

if(true)
    vfid    = get_vid_by_name(fileid,'fiducial');
    vcond1  = get_vid_by_name(fileid,'layer_conductivity');
    vcond2  = get_vid_by_name(fileid,'layer_conductivity_masked');
    vheight = get_vid_by_name(fileid,'height');
    
    %loop over lines
    figure(2);
    set(gcf,'defaultaxesfontsize',8);
    maximize_figure();    
    for k=1:5:nlines
        fid  = get_line(fileid,vfid,istart(k),icount(k));
        
        if(vcond1>=0)            
            cond1 = get_line(fileid,vcond1,[0 istart(k)],[30 icount(k)]);
        end
        
        if(vcond2>=0)
            cond2 = get_line(fileid,vcond2,[0 istart(k)],[30 icount(k)]);
        end
                
        if(vheight>=0)
            height = get_line(fileid,vheight,istart(k),icount(k));
        end
        
        clf;
        xl=[min(fid) max(fid)];        
        subplot(3,1,[1]);        
        if(vheight>0)
            plot(fid,height,'r');
        end
        grid on;
        xlim(xl);
        title(num2str(line(k)));
        
        subplot(3,1,[2:3]);
        cla;
        if(vcond1>0)
            semilogy(fid,cond1,'r');
        end
        hold on;
        
        if(vcond2>0)
            semilogy(fid,cond2,'b');
        end
        grid on;
        xlim(xl);                
        drawnow;
        pause(1)
    end
end
netcdf.close(fileid);
