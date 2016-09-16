
clc;
clear all;

%ncfile  = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/ncfiles/aem_conductivity/galeisbs.final.1e4.a.nc';
ncdir = 'Z:\projects\geophysics_netcdf\awags_levelled_conversion\ncfiles\';

%ncfile = [ncdir 'PIRSA_P830_TEiSA_M3_Mag.nc'];
%ncfile = [ncdir 'PIRSA_P833_SAEI_C4_tmi.nc'];
%ncfile = [ncdir 'PIRSA_P850_P851_TEiSA_A2_A3_B1_B2_Mag.nc'];
%ncfile = [ncdir 'PIRSA_P850_TEiSA_A1_Mag.nc'];
%ncfile = [ncdir 'GSNSW_P1027MAG.nc'];
%ncfile = [ncdir 'GSNSW_P754MAG.nc'];
%ncfile = [ncdir 'GSWA_P1235MAG.nc'];
%ncfile = [ncdir 'PIRSA_P853_TEiSA_D1_D2_D3_Mag_C.nc'];
%ncfile = [ncdir 'PIRSA_P853_TEiSA_D1_D2_D3_Mag_C.nc'];
%ncfile = [ncdir 'GSQP828MAG.nc'];
%ncfile = [ncdir 'P697MAG.nc']; 
%ncfile = [ncdir 'P452MAG.nc']; 
%ncfile = [ncdir 'GSQP1175MAG.nc']; 
ncfile = [ncdir 'P518MAG.nc']; 
%ncfile = [ncdir 'P406MAG.nc']; 

fileid = netcdf.open(ncfile,'NOWRITE');
istart = get_by_name(fileid,'_index_line');
icount = get_by_name(fileid,'_index_count');
x1  = get_by_name(fileid,'longitude_first');
x2  = get_by_name(fileid,'longitude_last');
y1  = get_by_name(fileid,'latitude_first');
y2  = get_by_name(fileid,'latitude_last');
lat = get_by_name(fileid,'latitude');
lon = get_by_name(fileid,'longitude');
tmi1 = get_by_name(fileid,'mag_awags');
bp  = get_by_name(fileid,'bounding_polygon');
line = get_by_name(fileid,'line');
nlines = length(istart);

bp(:,length(bp)+1)=bp(:,1);
%%
if(true)
    figure;
    maximize_figure();
    tic
    %plot(lon(1:1:end),lat(1:1:end),'.g');hold on;    
    plot([x1 x2]',[y1 y2]','b.'); hold on;
        
    %plot([x1 x2]',[y1 y2]','-b');
    plot(bp(1,:),bp(2,:),'-r.');
    %xlim([min([x1;x2]) max([x1;x2])]);    
    %ylim([min([y1;y2]) max([y1;y2])]);
    %xlim([min(lon) max(lon)]);
    %ylim([min(lat) max(lat)]);
    daspect([1 1 1]);
end

if(false);
    vfid   = get_vid_by_name(fileid,'fiducial');
    vtmi1  = get_vid_by_name(fileid,'mag_lev');
    vtmi2  = get_vid_by_name(fileid,'mag_awags');
    vheight = get_vid_by_name(fileid,'height');
    
    figure;
    set(gcf,'defaultaxesfontsize',8);
    maximize_figure();
    for k=1:100:nlines
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
        
        pause(0.2);
        pause;
    end
end
netcdf.close(fileid);
