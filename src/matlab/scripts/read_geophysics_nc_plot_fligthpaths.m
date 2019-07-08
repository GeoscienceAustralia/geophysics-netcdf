clc;
clear all;

plotactualflightpath = false;
plotboundingpolygon  = true;
plotstartendpoints   = false;

addpath('Z:\code\repos\geophysics-netcdf\src\matlab\functions\');

%basedir = 'Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\rb_working\';
basedir = 'Z:\projects\geophysics_netcdf\matlab_test\';
plotdir = [basedir];

%ncfiledir = 'Z:\projects\2017_AusAEM_Survey_01\final\line_data_em';
ncfiledir = 'Z:\projects\geophysics_netcdf\ncfiles\v2';
F = getfilelist(ncfiledir,'P441*.nc');
F = getfilelist(ncfiledir,'P1152*.nc');
%catalog  = 'catalog/uc0/rr2_dev/axi547/magnetic_line/catalog.html';
%F = get_ncfile_list(catalog);

k=1;
for i=1:1:1
    %ncfile = F(i).ncurl;
    ncfile = F(i).pathname;
    
    [p,n,e]  = fileparts(ncfile);
    filename = [n e];
    plotname = [plotdir filename '.jpg'];
    
    %if(exist(plotname,'file'))continue;end;
    disp([num2str(i) ' ' ncfile]);
    
    pt=0.1;    
    chunksize = 5e6;    
    [chosen_chunksize,fileid] = netcdf.open(ncfile,'NOWRITE',chunksize);
        
    [istart icount] = get_start_count(fileid);  
    nlines = length(istart);
    
              
    %% Figure
    close all;
    dark_figure(1);
    maximize_figure();  
    set(gcf,'InvertHardcopy','off');
    box on;
    hold on;
    
    if(plotactualflightpath)
        vx  = get_vid_by_stdname(fileid,'longitude');
        vy  = get_vid_by_stdname(fileid,'latitude');
        if(vx<0 || vy<0)
            vx  = get_vid_by_name(fileid,'Longitude');
            vy  = get_vid_by_name(fileid,'Latitude');
        end            
        if(vx<0 || vy<0)
            vx  = get_vid_by_name(fileid,'longitude');
            vy  = get_vid_by_name(fileid,'latitude');
        end            
                
        for k=1:1:nlines    
            %stride = floor(icount(k)/10);
            stride = icount(k)-1;
            x = get_line(fileid,vx,istart(k),icount(k),stride);
            y = get_line(fileid,vy,istart(k),icount(k),stride);
            plot(x,y,'-g+');
            daspect([1 1 1]);
            drawnow;
        end        
    end
    
    if(plotstartendpoints)
        x1 = get_by_name(fileid,'longitude_first');pause(pt);
        x2 = get_by_name(fileid,'longitude_last');pause(pt);
        y1 = get_by_name(fileid,'latitude_first');pause(pt);
        y2 = get_by_name(fileid,'latitude_last');pause(pt);               
        plot([x1 x2]',[y1 y2]','-b');        
    end
    
    if(plotboundingpolygon)
        bp = get_by_name(fileid,'bounding_polygon');                        
        if(~isempty(bp))
            patch(bp(1,:),bp(2,:),':','edgecolor','r','facecolor','none','linewidth',2);            
        end
    end
        
    daspect([1 1 1]);
    title(filename,'interpreter','none','fontsize',14);
    drawnow;
    saveas(gcf,plotname);     
    netcdf.close(fileid);
end

