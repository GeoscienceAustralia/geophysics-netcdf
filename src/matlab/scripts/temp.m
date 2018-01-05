ncfile  = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/AWAGS_Levelled_Line_Databases/mag_database_reformat_adjusted/netcdf/GSQ_P1248MAG.nc';
%ncfile  = 'Z:\projects\geophysics_netcdf\conversion_scripts\GSQ_P1248MAG.nc';

fileid = netcdf.open(ncfile,'NOWRITE');
istart = get_by_name(fileid,'index_line');
icount = get_by_name(fileid,'index_count');
line   = get_by_name(fileid,'line');
nlines  = length(istart);
vfid    = get_vid_by_name(fileid,'fiducial');
vx      = get_vid_by_name(fileid,'longitude');
vy      = get_vid_by_name(fileid,'latitude');
vtmi1   = get_vid_by_name(fileid,'mag_lev');
vtmi2   = get_vid_by_name(fileid,'mag_awags');
vheight = get_vid_by_name(fileid,'height');
for k=1:1:nlines
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
    
    %Remove an oder=porder polynomial trend - better padding extrapolation is required but this will do for now.
    porder = 2; pfit = polyfit(ind,mag,porder);
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
    cvd  = filter_coeff_vd(f,1);
    
    %Inverse FFT
    Sbp = cbp'.*S;%band pass filtered spectrem
    Svd = cvd'.*S;%1st vertical derivative filtered spectrem
    mag_bp = real(ifft(Sbp));%inverse bp spectrem
    mag_vd = real(ifft(Svd));%inverse vd spectrem    
end
netcdf.close(fileid);










