
clc;
clear all;

%ncfile  = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/ncfiles/aem_conductivity/galeisbs.final.1e4.a.nc';
ncdir = 'Z:\projects\geophysics_netcdf\awags_levelled_conversion\ncfiles\';
%ncfile = [ncdir 'PIRSA_P830_TEiSA_M3_Mag.nc'];
%ncfile = [ncdir 'PIRSA_P833_SAEI_C4_tmi.nc'];
ncfile = [ncdir 'PIRSA_P850_P851_TEiSA_A2_A3_B1_B2_Mag.nc'];
%ncfile = [ncdir 'PIRSA_P850_TEiSA_A1_Mag.nc'];

fileid = netcdf.open(ncfile,'NOWRITE');
istart = get_by_name(fileid,'_index_line');
icount = get_by_name(fileid,'_index_count');
x1  = get_by_name(fileid,'longitude_first');
x2  = get_by_name(fileid,'longitude_last');
y1  = get_by_name(fileid,'latitude_first');
y2  = get_by_name(fileid,'latitude_last');
lat = get_by_name(fileid,'latitude');
lon = get_by_name(fileid,'longitude');
bp  = get_by_name(fileid,'bounding_polygon');
line = get_by_name(fileid,'line');
nlines = length(istart);

%%
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
close all;

vfid   = get_vid_by_name(fileid,'fiducial');
vtmi1  = get_vid_by_name(fileid,'mag_mlev');
vtmi2  = get_vid_by_name(fileid,'mag_awags');
vheight = get_vid_by_name(fileid,'height');

figure;
set(gcf,'defaultaxesfontsize',8);
maximize_figure();
for k=1:100:nlines
    fid  = get_line(fileid,vfid,istart(k),icount(k));
    tmi1 = get_line(fileid,vtmi1,istart(k),icount(k));
    tmi2 = get_line(fileid,vtmi2,istart(k),icount(k));
    height = get_line(fileid,vheight,istart(k),icount(k));
    xl=[min(fid) max(fid)];
    
    if(true)  
        subplot(3,1,[1]);
        cla;
        plot(fid,height,'r');        
        grid on;
        xlim(xl);
        title(num2str(line(k)));
        
        subplot(3,1,[2:3]);
        cla;
        plot(fid,tmi1-mean_nonan(tmi1),'r');
        %plot(fid,tmi1,'r');
        hold on;
        plot(fid,100+tmi2-mean_nonan(tmi2),'b');
        grid on;
        xlim(xl);
    end
    pause(0.2);    
end
netcdf.close(fileid);
