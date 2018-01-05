
clc;
clear all;
addpath('Z:\code\repos\geophysics-netcdf\src\matlab');
basedir = 'Z:\projects\geophysics_netcdf\conversion_scripts\';
ncdir   = [basedir 'ncfiles\'];
ncfile = [ncdir 'GSSA_P1255MAG_Marree.nc'];
ncfile = [ncdir 'GSQ_P1248MAG.nc'];
ncfile = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_adjusted/netcdf/GSQ_P1248MAG.nc';
ncfile = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_2016_adjusted/netcdf/GSSA_P1255MAG_Marree.nc';

fileid = netcdf.open(ncfile,'NOWRITE');
pause(1);
istart = get_by_name(fileid,'index_line');
icount = get_by_name(fileid,'index_count');
line   = get_by_name(fileid,'line');
%lat = get_by_stdname(fileid,'latitude');
%lon = get_by_stdname(fileid,'longitude');

nlines = length(istart);

vfid    = get_vid_by_name(fileid,'fiducial');
vx   = get_vid_by_name(fileid,'longitude');
vy   = get_vid_by_name(fileid,'latitude');
vtmi1   = get_vid_by_name(fileid,'mag_lev');
vtmi2   = get_vid_by_name(fileid,'mag_awags');
vheight = get_vid_by_name(fileid,'height');

for k=1:floor(nlines/3):1
    fid  = get_line(fileid,vfid,istart(k),icount(k));
    
    x = get_line(fileid,vx,istart(k),icount(k));    
    y = get_line(fileid,vy,istart(k),icount(k));
    tmi_tielev = get_line(fileid,vtmi1,istart(k),icount(k));    
    tmi_awags  = get_line(fileid,vtmi2,istart(k),icount(k));    
    height = get_line(fileid,vheight,istart(k),icount(k));
    
    dx = 30 * 3600 * median(diff(x));
    dy = 30 * 3600 * median(diff(y));
    delta = sqrt(dx*dx+dy*dy);%average sample spacing in m
    
    
    %Choose data to work on
    mag = tmi_awags;
    N=2*floor(length(mag)/2);%even number
    mag    = tmi_awags(1:N);    
    ldist  = delta*(1:N);%line distance
    kmldist=ldist/1000;
    height = height(1:N);%height
    ind    = [1:N]';
           
    %Remove an oder=porder polynomial trend
    %better padding extrapolation is required but this will do for now.
    porder = 2;    
    pfit = polyfit(ind,mag,porder);
    mag_detrend = mag - polyval(pfit,ind);        
        
    %FFT data size and frequencies    
    v = mag_detrend;
    S = fft(v);%the spectrem    
    f = fft_frequency(length(S+1),delta);       
        
    %Determine band-pass filter coefficients             
    lcut_wl = 10000;%low cut wavelength in m (note lcut_wl>hcut_wl)
    hcut_wl = 1000;%hcut cut wavelength in m
    lcut_f = 1./lcut_wl;%low cut frequency in cycles/m
    hcut_f = 1./hcut_wl;%high cut frequency in cycles/m    
    cbp = filter_coeff_bp_bw(f,lcut_f,hcut_f,1);    
    %chp = filter_coeff_hp_bw(f,lcut_f,1);
    %clp = filter_coeff_lp_bw(f,hcut_f,1);
    cvd  = filter_coeff_vd(f,1);
    
    %Inverse FFT
    Sbp = cbp'.*S;%band pass filtered spectrem
    Svd = cvd'.*S;%1st vertical derivative filtered spectrem
    mag_bp = real(ifft(Sbp));%inverse bp spectrem
    mag_vd = real(ifft(Svd));%inverse vd spectrem
        
    %%
    figure(1);  
    clf;
    maximize_figure();
    h1=loglog(f,abs(S),'b');    
    hold on;
    h2=loglog(f,abs(Sbp),'r');        
    h3=loglog(f,abs(Svd),'b');        
    grid on; box on;
    ylabel('Power');    
    
    %xlabel('Frequency (cycles/m)');        
    %Plot as wavelength rather than frequency    
    xt=get(gca,'xtick');
    set(gca,'xticklabel',1./xt);
    xlabel('Wavelength (m)'); 
    lh=legend([h1 h2 h3],'detrended', 'band pass','vd');
    drawnow;
    pause(0.1);
    
    %%    
    figure(2);
    clf;
    maximize_figure();
    set(gcf,'defaultaxesfontsize',8);        
    xlimits=[min(kmldist) max(kmldist)];    
    
    %Top panel height    
    subplot(5,1,[1]);    
    plot(kmldist,height,'k');    
    grid on;
    xlim(xlimits);
    title(num2str(line(k)));
    ylabel('Height (m)');
    
    %Bottom panel magnetics
    subplot(5,1,[2:3]);        
    h1=plot(kmldist,mag_detrend,'k');
    hold on;
    h2=plot(kmldist,mag_bp,'r');
    grid on;
    xlim(xlimits);    
    xlabel('Distance (km)');
    ylabel('TMI (nT)');
    lh=legend([h1 h2],'detrended', 'band pass');
    set(lh,'location','northeast');    
    
    %Bottom panel magnetics
    subplot(5,1,[4:5]);            
    h1=plot(kmldist,mag_vd,'b');    
    grid on;
    xlim(xlimits);
    ylim([-0.01 0.01]);
    xlabel('Distance (km)');
    ylabel('TMI (nT)');    
    lh=legend([h1],'vd');
    set(lh,'location','northeast');    
    drawnow;
    pause(1);
end
netcdf.close(fileid);



