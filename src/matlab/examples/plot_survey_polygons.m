clc;
clear all;

addpath('..\');

basedir = 'Z:\projects\geophysics_netcdf\conversion_scripts\';
ncfiledir = 'Z:\projects\geophysics_netcdf\ncfiles\v2\';
F = getfilelist(ncfiledir,'P441*.nc');
%F = getfilelist(ncfiledir,'P1152*.nc');
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
BP=[];
for i=1:1:length(F)
    %ncfilepath = F(i).ncurl;            
    ncfilepath  = F(i).pathname;    
    N = GeophysicsNCFile(ncfilepath);                
    bp = N.get_bounding_polygon();                            
    if(isempty(bp))              
       msg = sprintf('%s has no bounding polygon',N.path);
       warning(msg);       
    else
       BP(k).poly = bp;
       k=k+1;                
    end
    N.close();
end

%%

dark_figure(1);

maximize_figure();
for i=1:1:length(BP);
    p = BP(i).poly;
    patch(p(1,:),p(2,:),':','edgecolor','r','facecolor','none','linewidth',2);                
    hold on; grid on; box on; drawnow;
end
daspect([1 1 1]);
%plotname = [plotdir '000_overall.jpg'];
%saveas(gcf,plotname);


