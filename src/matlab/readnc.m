clc;
clear all;

%ncfile   = '..\ncfiles\magnetics\GSNSW_P500MAG.nc';
%ncfile   = '..\ncfiles\magnetics\GSQP1029MAG.nc';
ncfile   = 'http://dapds00.nci.org.au/thredds/dodsC/uc0/rr2_dev/rcb547/magrad_tests_indexed_v2/GSQP1029MAG.nc';

tic
ncid = netcdf.open(ncfile,'NOWRITE');
disp(ncid)
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);
disp('Inquire');
toc

tic
for vi=0:1:numvars-1    
    [varname, xtype, dimids, numatts] = netcdf.inqVar(ncid,vi);    
    for ai=0:1:numatts-1         
        attname = netcdf.inqAttName(ncid,vi,ai);
        attval  = netcdf.getAtt(ncid,vi,attname);                   
        if(strcmp(attname,'my_standard_name'))            
            if(strcmp(attval,'fiducial_number'))
                vid_fid = vi;            
            elseif(strcmp(attval,'index_lines'))
                vid_index_lines = vi;            
            elseif(strcmp(attval,'longitude_first'))
                vid_x1 = vi;            
            elseif(strcmp(attval,'longitude_last'))
                vid_x2 = vi;            
            elseif(strcmp(attval,'latitude_first'))
                vid_y1 = vi;            
            elseif(strcmp(attval,'latitude_last'))
                vid_y2 = vi;            
            elseif(strcmp(attval,'total_magnetic_intensity'))
                vid_tmi1 = vi;            
            elseif(strcmp(attval,'total_magnetic_intensity_tielevelled'))
                vid_tmi2 = vi;            
            elseif(strcmp(attval,'line_number'))
                vid_line = vi;
            end            
        end        
    end
end
disp('Var ids');
toc

tic
index  = netcdf.getVar(ncid,vid_index_lines);
nlines = length(index)-1;
start  = index(1:nlines);
count  = diff(index);
line   = netcdf.getVar(ncid,vid_line);
disp('Indexes');
toc

tic
x1 = netcdf.getVar(ncid,vid_x1);
x2 = netcdf.getVar(ncid,vid_x2);
y1 = netcdf.getVar(ncid,vid_y1);
y2 = netcdf.getVar(ncid,vid_y2);
disp('End points');
toc

%%
figure;
maximize_figure();
tic
plot([x1 x2]',[y1 y2]','-b');
disp('Render lines');
toc
daspect([1 1 1]);

tic
for k=1:1000:nlines    
     fid  = netcdf.getVar(ncid,vid_fid,start(k),count(k));
     tmi1 = netcdf.getVar(ncid,vid_tmi1,start(k),count(k));
     tmi2 = netcdf.getVar(ncid,vid_tmi2,start(k),count(k));
     pause(0.01);
end
disp('Get 3 variables');
toc

%%
figure;
maximize_figure();
for k=1:1000:nlines
    fid  = netcdf.getVar(ncid,vid_fid,start(k),count(k));
    tmi1 = netcdf.getVar(ncid,vid_tmi1,start(k),count(k));
    tmi2 = netcdf.getVar(ncid,vid_tmi2,start(k),count(k));
    
    if(true)
        tmi1(tmi1==-Inf)=NaN;
        tmi2(tmi2==-Inf)=NaN;

        cla;
        plot(fid,tmi1-mean_nonan(tmi1),'r');
        hold on;
        plot(fid,tmi2-mean_nonan(tmi2),'b');    
        title(num2str(line(k)));        
    end
    pause(0.01);    
end
disp('Get 3 variables');
toc
netcdf.close(ncid);





