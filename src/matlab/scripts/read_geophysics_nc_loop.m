
clc;
clear all;

addpath('Z:\code\repos\geophysics-netcdf\src\matlab');

catalog1 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/rad_database_reformat_adjusted/netcdf/catalog.html'; 
catalog2 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/rad_database_reformat_2016_adjusted/netcdf/catalog.html'; 

%basedir = 'Z:\projects\geophysics_netcdf\conversion_scripts\';
%ncdir   = [basedir 'ncfiles\'];
%plotdir = [basedir 'plots\'];

k=1;
for si=1:1:1    
    ncfile = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_adjusted/netcdf/GSQ_P1248MAG.nc';
                     
    fileid = netcdf.open(ncfile,'NOWRITE');    
        
    istart = get_by_name(fileid,'index_line');
    icount = get_by_name(fileid,'index_count');
    x1  = get_by_name(fileid,'longitude_first');
    x2  = get_by_name(fileid,'longitude_last');
    y1  = get_by_name(fileid,'latitude_first');
    y2  = get_by_name(fileid,'latitude_last');
    %lat = get_by_stdname(fileid,'latitude');
    %lon = get_by_stdname(fileid,'longitude');
    bp     = get_by_name(fileid,'bounding_polygon');
    line   = get_by_name(fileid,'line');
    nlines = length(istart);
    
    bp(:,length(bp)+1)=bp(:,1);
    
    BP(k).poly = bp;
    k=k+1;
    
    if(false)
        close all;
        figure;
        maximize_figure();
        tic
        plot([x1 x2]',[y1 y2]','-b');
        hold on;
        
        %plot(lon(1:100:end),lat(1:100:end),'g.');
        plot([x1 x2]',[y1 y2]','-b');
        plot(bp(1,:),bp(2,:),'-r.');
        %xlim([min([x1;x2]) max([x1;x2])]);
        %ylim([min([y1;y2]) max([y1;y2])]);
        daspect([1 1 1]);
        
        [p,n,e] = fileparts(ncfile);
        filename = [n e];
        title(filename,'interpreter','none','fontsize',8);                
        plotname = [plotdir filename '.jpg'];
        saveas(gcf,plotname);
    end
        
    if(true);
        vfid    = get_vid_by_name(fileid,'fiducial');
        vtmi1   = get_vid_by_name(fileid,'mag_lev');
        vtmi2   = get_vid_by_name(fileid,'mag_awags');
        vheight = get_vid_by_name(fileid,'height');
        
        figure;
        set(gcf,'defaultaxesfontsize',8);
        maximize_figure();
        for k=1:10:nlines
            fid  = get_line(fileid,vfid,istart(k),icount(k));
            
            if(vtmi1>=0)
                tmi1 = get_line(fileid,vtmi1,istart(k),icount(k));
            end
            if(vtmi2>=0)
                tmi2 = get_line(fileid,vtmi2,istart(k),icount(k));
            end
            
            if(vheight>=0)
                height = get_line(fileid,vheight,istart(k),icount(k));
            end
            xl=[min(fid) max(fid)];
            
            subplot(3,1,[1]);
            cla;
            if(vheight>0)
                plot(fid,height,'r');
            end
            grid on;
            xlim(xl);
            title(num2str(line(k)));
            
            subplot(3,1,[2:3]);
            cla;
            if(vtmi1>0)
                plot(fid,tmi1-mean_nonan(tmi1),'r');
            end
            hold on;
            
            if(vtmi2>0)
                plot(fid,100+tmi2-mean_nonan(tmi2),'b');
            end
            grid on;
            xlim(xl);                        
            pause(0.001);            
        end
    end        
    netcdf.close(fileid);        
end
