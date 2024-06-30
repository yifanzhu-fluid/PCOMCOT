function plotBathymetry(layerID)

fn = sprintf('PCOMCOToutput/_bathymetry%02d.dat',layerID);
[x, y, b] = COMCOT_readBinaryDataSnapshot(fn);
dx = x(2)-x(1); dy = y(2)-y(1);

b(b(:,:)<=0.0) = nan;


close all
figure; hold on
set(gcf,'Units','normal','Position',[0.1,0.1,0.8,0.8]);

im = imagesc(x,y,b);
set(im,'alphadata',~isnan(b));
colormap(flipud(jet))
colorbar;

if(exist('Stations.ctl','file') == 2)
    [StaInfo, nStations] = readStationsInfoCTL('Stations.ctl');
    for i = 1:nStations
        x1 = StaInfo(i).lon; y1 = StaInfo(i).lat;
        plot(x1,y1,'Color','k','Marker','^','MarkerSize',8,'MarkerFaceColor','k');
        text(x1-15*dx,y1-15*dy,StaInfo(i).name);
    end
end

ilayer = 1;
while 1>0
    ilayer = ilayer+1; fn = sprintf('_bathymetry%02d.dat',ilayer);
    if(exist(fn)~=2); break; end
    xtmp = load(['_xcoordinate',sprintf('%02d',ilayer),'.dat']);
    ytmp = load(['_ycoordinate',sprintf('%02d',ilayer),'.dat']);
    x1 = min(xtmp); x2 = max(xtmp); y1 = min(ytmp); y2 = max(ytmp);
    plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'-w','LineWidth',2);
    text(x1+10*dx,y2-20*dy,sprintf('layer%02d',ilayer),'Color','w','FontSize',15,'FontWeight','bold');
end
set(gca,'FontSize',15); axis equal; 
ylim([min(y),max(y)]);
xlim([min(x),max(x)]);
