function plotSnapshot(fn)

if(nargin == 0)
    fn = 'z_01_000000.dat';
end

[x,y,b] = COMCOT_readBinaryDataSnapshot(fn);
i = strfind(fn,filesep);
if isempty(i)
    dir = ''; dataFile = fn;
else
    dir = fn(1:i(end)); dataFile = fn(i(end)+1:end);
end
i = strfind(dataFile,'_');
ilayer = dataFile(i(1)+1:i(1)+2);
[x,y,h] = COMCOT_readBinaryDataSnapshot([dir,'_bathymetry',ilayer,'.dat']);
if ~(contains(dataFile,'M')||contains(dataFile,'N'))
    b(h<=0&b+h<=0) = nan;
end

close all
figure; hold on
set(gcf,'Units','normal','Position',[0.3,0.3,0.6,0.6]);

bmax = max(max(b)); bmin = min(min(b));
im = imagesc(x,y,b,[bmin,bmax]);
contour(x,y,h,[0,0],'k-','LineWidth',1.0);
set(gca,'YDir','normal');
colormap('jet');
set(im,'alphadata',~isnan(b));
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],'k-','LineWidth',1.5);
colorbar;
title(dataFile,'interpreter','none');
axis equal;
xlim([min(x),max(x)]);
ylim([min(y),max(y)]);




