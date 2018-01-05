clc;
clear all;

plotactualflightpath = true;
plotstartendpoints   = true;
plotboundingpolygon  = true;

addpath('Z:\code\repos\geophysics-netcdf\src\matlab\functions\');

basedir = 'Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\rb_working\';
plotdir = [basedir 'awags_conversions\plots\awags_plots\'];

ncfiledir = 'Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\awags_survey_reformat\netcdf\';
F = getfilelist(ncfiledir,'*.nc');

%catalog1 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_adjusted/netcdf/catalog.html';
%catalog2 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/catalog.html';
%plotdir = [basedir 'mag_plots\'];

%catalog1 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/rad_database_reformat_adjusted/netcdf/catalog.html';
%catalog2 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/rad_database_reformat_2016_adjusted/netcdf/catalog.html';
%plotdir = [basedir 'rad_plots\'];
%F1 = get_ncfile_list(catalog1);
%F2 = get_ncfile_list(catalog2);
%F = [F1 F2];

%catalog1 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/awags_survey_reformat/netcdf/catalog.html';
%F = get_ncfile_list(catalog1);

k=1;
for i=1:1:length(F)
    %ncfile = F(i).ncurl;
    ncfile = F(i).pathname;
    
    [p,n,e] = fileparts(ncfile);
    filename = [n e];
    plotname = [plotdir filename '.jpg'];
    
    %if(exist(plotname,'file'))continue;end;
    disp([num2str(i) ' ' ncfile]);
    
    pt=0.1;
    chunksize = 20e6;
    [chosen_chunksize,fileid] = netcdf.open(ncfile,'NOWRITE',chunksize);
    pause(pt);
    
    istart = get_by_name(fileid,'index_line');pause(pt);
    icount = get_by_name(fileid,'index_count');pause(pt);
    line   = get_by_name(fileid,'line');pause(pt);
    nlines = length(istart);pause(pt);
    
    %figure
    close all;
    figure;
    maximize_figure();  
    box on;
    hold on;
    
    if(plotactualflightpath)
        vx  = get_vid_by_name(fileid,'longitude');
        vy  = get_vid_by_name(fileid,'latitude');
        for k=1:1:nlines
            x = get_line(fileid,vx,istart(k),icount(k));
            y = get_line(fileid,vy,istart(k),icount(k));
            plot(x(1:100:end),y(1:100:end),'g.');
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
        plot(bp(1,:),bp(2,:),'-r.');        
    end
        
    daspect([1 1 1]);
    title(filename,'interpreter','none','fontsize',8);
    saveas(gcf,plotname); 
    
    netcdf.close(fileid);
end

