clc;
clear all;

addpath('Z:\code\repos\geophysics-netcdf\src\matlab');
basedir = 'Z:\projects\geophysics_netcdf\conversion_scripts\';

ncfiledir = 'Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\awags_survey_reformat\netcdf\';
F = getfilelist(ncfiledir,'*.nc');

%catalog1 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_adjusted/netcdf/catalog.html'; 
%catalog2 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/catalog.html'; 
%plotdir = [basedir 'mag_plots\'];

%catalog1 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/rad_database_reformat_adjusted/netcdf/catalog.html'; 
%catalog2 = 'catalog/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/rad_database_reformat_2016_adjusted/netcdf/catalog.html'; 
%plotdir  = [basedir 'rad_plots\'];

%F1 = get_ncfile_list(catalog1);
%F2 = get_ncfile_list(catalog2);
%F = [F1 F2];

k=1;
for i=1:1:length(F)
    %ncfile = F(i).ncurl;            
    ncfile = F(i).pathname;
    [p,n,e] = fileparts(ncfile);
    filename = [n e];            
    
    disp([num2str(i) ' ' ncfile]);
    
    fileid = netcdf.open(ncfile,'NOWRITE');        
    istart = get_by_name(fileid,'index_line');
    icount = get_by_name(fileid,'index_count');
        
    bp     = get_by_name(fileid,'bounding_polygon');
    line   = get_by_name(fileid,'line');
    nlines = length(istart);
    
    bp(:,length(bp)+1)=bp(:,1);
    BP(k).poly = bp;
    k=k+1;        
    netcdf.close(fileid);
    pause(0.5);    
end

%%
figure;
maximize_figure();
for i=1:1:length(BP);
    p = BP(i).poly;
    plot(p(1,:),p(2,:));
    hold on;
    grid on;
    box on;
    drawnow;    
end
daspect([1 1 1]);
%plotname = [plotdir '000_overall.jpg'];
%saveas(gcf,plotname);


