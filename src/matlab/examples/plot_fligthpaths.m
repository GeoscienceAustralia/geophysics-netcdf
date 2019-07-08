clc;
clear all;

plotactualflightpath = true;
plotboundingpolygon  = true;
plotstartendpoints   = true;

addpath('..\');

%basedir = 'Y:\ops\gap\geophysical_methods\mag_rad\AWAGS_Levelled_Databases\rb_working\';
basedir = 'Z:\projects\geophysics_netcdf\matlab_test\';
plotdir = [basedir];

%ncfiledir = 'Z:\projects\2017_AusAEM_Survey_01\final\line_data_em';
ncfiledir = 'Z:\projects\geophysics_netcdf\ncfiles\v2';
F = getfilelist(ncfiledir,'P441*.nc');
%F = getfilelist(ncfiledir,'P1152*.nc');
%catalog  = 'catalog/uc0/rr2_dev/axi547/magnetic_line/catalog.html';
%F = get_ncfile_list(catalog);

k=1;
for i=1:1:1
    %ncfilepath = F(i).ncurl;
    ncfilepath  = F(i).pathname;
    
    [p,n,e]  = fileparts(ncfilepath);
    filename = [n e];
    plotname = [plotdir filename '.jpg'];
    
    %if(exist(plotname,'file'))continue;end;
    disp([num2str(i) ' ' ncfilepath]);
                
    N = GeophysicsNCFile(ncfilepath);
    %N.init_start_count();
    
    nlines = N.nlines();
    var1 = N.get_var_by_name('longitude');        
    var2 = N.get_var_by_attval('standard_name','longitude');        
    var3 = N.get_var_by_stdname('longitude');
    var4 = N.get_var_by_longname('total_magnetic_intensity_anomaly_tie_line_levelled');    
        
    %% Figure
    close all;
    dark_figure(1);
    maximize_figure();  
    set(gcf,'InvertHardcopy','off');
    box on;
    hold on;
    
    if(plotactualflightpath)
        vx  = N.get_var_by_stdname('longitude');
        vy  = N.get_var_by_stdname('latitude');
        
        if(isempty(vx) || isempty(vy))
            vx  = N.get_var_by_name('longitude');
            vy  = N.get_var_by_name('latitude');
        end            
        
        if(isempty(vx) || isempty(vy))
            vx  = N.get_var_by_name('Longitude');
            vy  = N.get_var_by_name('Latitude');
        end            
                        
        for k=1:1:nlines    
            stride = floor(N.nsamples(k)/10);
            %stride = N.nsamples(k) - 1;
            stride = 1;
            x = N.get_vardata_lineindex(vx,k,stride);
            y = N.get_vardata_lineindex(vy,k,stride);
            plot(x,y,'-g');
            daspect([1 1 1]);
            drawnow;
        end        
    end
    
    if(plotstartendpoints)
        [xstart,xend,ystart,yend] = N.get_line_start_end_points();
        plot([xstart xend]',[ystart yend]','-b');
    end
    
    if(plotboundingpolygon)
        bp = N.get_bounding_polygon();                        
        if(~isempty(bp))
            patch(bp(1,:),bp(2,:),':','edgecolor','r','facecolor','none','linewidth',2);            
        end
    end
    
    N.close();
    
    daspect([1 1 1]);
    title(filename,'interpreter','none','fontsize',14);
    drawnow;
    saveas(gcf,plotname);     
    
end

