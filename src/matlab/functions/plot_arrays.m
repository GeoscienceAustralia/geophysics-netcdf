function plotme()
    clc;
    clear all;
    close all;
    
    red   = [1.0 0.5 0.5];    
    green = [0.5 1.0 0.5];
    blue  = [0.5 0.5 1.0];
    yellow = [0.5 1.0 1.0];
    
    fw=27.5;    
    fh=18;
    
    
    default_figure(fw,fh,12);    
    set(gcf,'position',[2.1167 0 fw fh]);
    set(gcf,'defaultaxeslinewidth',1);
    set(gcf,'defaultlinelinewidth',1);
    set(gcf,'menubar','none');
    set(gcf,'toolbar','none');
    
    
    %pos = [2.1167 0 27.5167 16.7746];
    pos = [0 0 fw fh];
    axes('position',pos,'units','centimeters');        
    
    set(gca,'layer','top');
    xlabel([]); ylabel([]);
    set(gca,'xtick',[]); set(gca,'ytick',[]);                
    set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);                
    hold on;
    
    daspect([1 1 1]);    
    xlim([0 fw]);
    ylim([0 fh]);
    
    ns = 80;
    nl = 10;
    nlay = 4;
    nwin = 12;
    nwave=7;
    pv=18;
    po=2;
    y=fh-0.2;
    dy=0.6;
    x=4;
    
    plotarray('point' ,x,y,1,ns,red);y=y-dy;
    plotarray('line',x,y,1,nl,red);y=y-dy;
    plotarray('waveform',x,y,1,nwave,red);y=y-dy;
    plotarray('window',x,y,1,nwin,red);y=y-dy;
    plotarray('layer',x,y,1,nlay,red);y=y-dy;
    plotarray('polygonvertex',x,y,1,pv,red);y=y-dy;
    plotarray('polygonordinate',x,y,1,po,red);y=y-dy;
    
        
    x= 13;
    y= fh-1;
    plotarray('project',x,y,1,nl,blue);y=y-dy;    
    plotarray('flight',x,y,1,nl,blue);y=y-dy;    
    plotarray('date',x,y,1,nl,blue);y=y-dy;   
    plotarray('bearing',x,y,1,nl,blue);y=y-dy;
    plotarray('longitude_first',x,y,1,nl,blue);y=y-dy;
    plotarray('longitude_last',x,y,1,nl,blue);y=y-dy;
    plotarray('latitude_first',x,y,1,nl,blue);y=y-dy;
    plotarray('latitude_last',x,y,1,nl,blue);y=y-dy;
    plotarray('index_count',x,y,1,nl,yellow);y=y-dy;    
    xind=x; yind=y;
    plotarray('index_start',x,y,1,nl,yellow);y=y-dy-2;
    
    
    x= 21;
    y= fh-1;
    clr = [0.8 0.8 0.8]; 
    plotarray('window_start_time',x,y,1,nwin,clr); y=y-dy;
    plotarray('window_end_time',x,y,1,nwin,clr); y=y-dy;
    plotarray('waveform_time',x,y,1,nwave,clr); y=y-dy;
    plotarray('waveform_current',x,y,1,nwave,clr); y=y-dy;    
    plotarray('bounding_polygon',x,y,po,pv,clr);y=y-dy;
    
    
    x = 4;
    y = 10;
    plotarray('fiducial',x,y,1,ns,green);y=y-dy;        
    plotarray('latitude' ,x,y,1,ns,green);y=y-dy;
    plotarray('longitude' ,x,y,1,ns,green);y=y-dy;
    plotarray('magnetics' ,x,y,1,ns,green);y=y-dy;    
    plotarray('EM_Z',x,y,nwin,ns,green);y=y-dy-(nwin*0.3);    
    plotarray('conductivity',x,y,nlay,ns,green);y=y-dy-(nlay*0.3);
    plotarray('thickness' ,x,y,nlay,ns,green);y=y-dy-(nlay*0.3);
    
    ds=[0 3 7 10 20 30 40 46 65 70]
    xdst = 4+0.3/2;
    ydst = 10;
    yind = yind-0.3/2;
    xind = xind-0.3/2;
    for i=1:1:nl
        px = [xind+i*0.3 xdst+0.3*ds(i)]/fw;
        py = [yind       ydst]/fh;
        h=annotation('arrow',px,py);    
        set(h,'linewidth',1);
        set(h,'headlength',3);
        set(h,'headwidth',3);
    end
    
end

function plotarray(label,x0,y0,m,n,clr)
    dx=.3; dy=.3;
    for i=1:1:m
        y = y0 - (i-1)*dy;
        for j=1:1:n
            x = x0 + (j-1)*dx;
            %plot([x x x+dx x+dx x],[y y-dy y-dy y y],'k');            
            hold on;
            p=patch([x x x+dx x+dx],[y y-dy y-dy y],'k');            
            set(p,'facecolor',clr);
            set(p,'edgecolor',[0 0 0]);
        end        
    end    
    tx = x0 - 0.2;
    ty = y0 - (m/2)*dy;
    th = text(tx,ty,label);    
    set(th,'units','centimeters');
    set(th,'verticalalignment','middle');
    set(th,'horizontalalignment','right');
    set(th,'interpreter','none');
end




