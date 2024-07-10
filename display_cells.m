function display_cells(metadata,stack,CELLS,LBL_FijiROI)
%DISPLAY_CELLS Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 27-Mar-2023 NL

%% 

close all
hfig = figure;
ax = axes('Parent',hfig);

im = mean(squeeze(stack.Ch2),3);
im = uint16(im);
im = imadjust(im,stretchlim(im,[0.01 0.9999]));
im_d = imdilate(LBL_FijiROI,strel("square",3));
im_e = imerode(LBL_FijiROI,strel("square",3));
imROI = imsubtract(im_d,im_e);
imshow(labeloverlay(im,imROI))
% imshow(LBL_FijiROI,[])
ax.Title.String = metadata.FileName;
ax.Title.Interpreter = 'none';
ax.Subtitle.String = 'tdTomato time average';
ax.Visible = 'on';
ax.XTick = [];
ax.YTick = [];
hold(ax,'on')
t = metadata.time;
for idxC = 1:height(CELLS)
    xy = CELLS.Centroid_xy_pix{idxC};
    x = nan(metadata.DimT,1);
    y = nan(metadata.DimT,1);
    x(CELLS.TimeFrame{idxC}) = xy(:,1);
    y(CELLS.TimeFrame{idxC}) = xy(:,2);
    x = [x;nan];
    y = [y;nan];
    [r,c] = ind2sub([metadata.DimY, metadata.DimX],find(LBL_FijiROI==CELLS.TrackID(idxC),1,'First'));
%     text(xy(end,1),xy(end,2),num2str(CELLS.TrackID(idxC),'%02.0f'),...
%         'Color',[1 1 1],'BackgroundColor',[0.2 0.2 0.2])
    text(c,r,num2str(CELLS.TrackID(idxC),'%02.0f'),...
        'Color',[1 1 1],'BackgroundColor',[0.2 0.2 0.2],...
        'FontSize',10)
    hp = patch(x,y,seconds([t;nan]),'EdgeColor','interp','LineWidth',1);  

end
hc = colorbar;
hc.Label.String = 'time (s)';
hold(ax,'off')
hfig.WindowState = 'maximized';
drawnow
exportgraphics(ax,fullfile(metadata.FilePath,['Cells kept with Fiji.png']),'ContentType','auto')
% exportgraphics(ax,fullfile(metadata.FilePath,[metadata.FileName '.png']),'ContentType','auto')

%%

for idxC = 1:height(CELLS)
vid_obj = VideoWriter(fullfile(metadata.FilePath,['Cell ID ' num2str(CELLS.TrackID(idxC)) '.mp4']),'MPEG-4');
vid_obj.Quality = 100;
vid_obj.FrameRate = 1/seconds(metadata.ResT);
open(vid_obj);


close all
hfig = figure('WindowState','maximized');
htl  = tiledlayout(6,3,"TileSpacing","tight","Padding","tight",'Parent',hfig);
htl.Title.String = [metadata.FileName ' cell: ' num2str(CELLS.TrackID(idxC),'%02.0f')];
htl.Title.Interpreter = 'none';

ax(1) = nexttile(htl,1,[1 1]);
ax(2) = nexttile(htl,4,[1 1]);
ax(3) = nexttile(htl,7,[1 1]);
ax(4) = nexttile(htl,10,[1 1]);
ax(5) = nexttile(htl,13,[1 1]);
ax(6) = nexttile(htl,16,[1 1]);

ax(7) = nexttile(htl,2,[3 1]);
ax(8) = nexttile(htl,3,[3 1]);
ax(9) = nexttile(htl,11,[3 1]);
ax(10) = nexttile(htl,12,[3 1]);

t = metadata.time;

y = nan(metadata.DimT,1);
y(CELLS.TimeFrame{idxC}) = CELLS.Salsa6fSig{idxC};
hl_salsa = line(t,y,'Parent',ax(1));
ax(1).YLabel.String = 'Salsa6f';
% ax(1).XTick = [];
ax(1).XGrid = 'on';
ax(1).Box = 'on';

y = nan(metadata.DimT,1);
y(CELLS.TimeFrame{idxC}) = CELLS.Area{idxC}*metadata.ResX*metadata.ResY;
hl_area = line(t,y   ,'Parent',ax(2));
ax(2).YLabel.String = 'Area (\mum^2)';
% ax(2).XTick = [];
ax(2).XGrid = 'on';
ax(1).Box = 'on';

y = nan(metadata.DimT,1);
y(CELLS.TimeFrame{idxC}) = CELLS.Circularity{idxC};
hl_circ = line(t,y,'Color',[0 0.4470 0.7410],...
    'DisplayName','Circularity', 'Parent',ax(3));
% y = nan(metadata.DimT,1);
% y(CELLS.TimeFrame{idxC}) = CELLS.Solidity{idxC};
% line(t,y,'Color',[0.8500 0.3250 0.0980],...
%     'DisplayName','Solidity','Parent',ax(3));
ax(3).YLabel.String = 'Circularity (-)';
% ax(3).XTick = [];
ax(3).XGrid = 'on';
% legend(ax(3),'Location','southeast')
ax(1).Box = 'on';


y = nan(metadata.DimT,1);
y(CELLS.TimeFrame{idxC}) = CELLS.ShapeRatio{idxC};
hl_shrt = line(t,y,'Color',[0 0.4470 0.7410],...
   'DisplayName','Shape ratio','Parent',ax(4));
% y = nan(metadata.DimT,1);
% y(CELLS.TimeFrame{idxC}) = CELLS.Eccentricity{idxC};
% line(t,y,'Color',[0.8500 0.3250 0.0980],...
%     'DisplayName','Eccentricity','Parent',ax(4));
ax(4).YLabel.String = 'Shape ratio (-)';
% ax(4).XTick = [];
ax(4).XGrid = 'on';
% legend(ax(4),'Location','southeast')
ax(4).Box = 'on';

% y = nan(metadata.DimT,1);
y = CELLS.speed{idxC};
hl_speed = line(t,y,'Parent',ax(5));
ax(5).YLabel.String = 'Speed (\mum\cdots^{-1})';
ax(5).XGrid = 'on';
ax(5).Box = 'on';

Tclip = 0.25;
x = CELLS.lag{idxC}(2:round(length(CELLS.lag{idxC})*Tclip));
y = CELLS.MSD{idxC}(2:round(length(CELLS.lag{idxC})*Tclip));
hl_msd = line(x,y,'Parent',ax(6));       
ax(6).XLabel.String = 'lag \tau(s)';
ax(6).YLabel.String = 'MSD (\mum^2)';
ax(6).XGrid = 'on';
ax(6).Box = 'on';
ax(6).Title.String = ['MSD = D\tau^\alpha+C - ',...
      'D = ' num2str(CELLS.MSD_fit{idxC}.D,'%2.2e') '\mum^2/s, \alpha = ' num2str(CELLS.MSD_fit{idxC}.alpha,'%2.1f')];


idxT = 1;

h_rp(1) = line(hl_salsa.XData(idxT),hl_salsa.YData(idxT),...
    'Marker','o','Parent',ax(1));
h_rp(2) = line(hl_area.XData(idxT),hl_area.YData(idxT),...
    'Marker','o','Parent',ax(2));
h_rp(3) = line(hl_circ.XData(idxT),hl_circ.YData(idxT),...
    'Marker','o','Parent',ax(3));
h_rp(4) = line(hl_shrt.XData(idxT),hl_shrt.YData(idxT),...
    'Marker','o','Parent',ax(4));
h_rp(5) = line(hl_speed.XData(idxT),hl_speed.YData(idxT),...
    'Marker','o','Parent',ax(5));


bw = false(metadata.DimY,metadata.DimX,metadata.DimT);
PixelIdxList = cat(1,CELLS.PixelIdxList{idxC}{:});
bw(PixelIdxList) = true;
bw = any(bw,3);
[row,col] = ind2sub([metadata.DimY metadata.DimX],find(bw));
y_lim     = [min(row) max(row)];
SZheight  = y_lim(2)-y_lim(1);
x_lim     = [min(col) max(col)];
SZwidth   = x_lim(2)-x_lim(1);
extbox = 0.2;
y_lim = [max(round(y_lim(1)-extbox/2*SZheight),1):min(round(y_lim(2)+extbox/2*SZheight),metadata.DimY)];
x_lim = [max(round(x_lim(1)-extbox/2*SZwidth) ,1):min(round(x_lim(2)+extbox/2*SZwidth) ,metadata.DimX)];

Ch1 = squeeze(stack.Ch1(y_lim,x_lim,1,:));
Ch2 = squeeze(stack.Ch2(y_lim,x_lim,1,:));
him(1) = imshow(Ch1(:,:,idxT),[min(Ch1(:)) max(Ch1(:))],...
    'Parent',ax(7));
him(2) = imshow(Ch2(:,:,idxT),[min(Ch2(:)) max(Ch2(:))],...
    'Parent',ax(8));
Ch3 = nan(metadata.DimY,metadata.DimY,metadata.DimT);
Ch3(cat(1,CELLS.PixelIdxList{idxC}{:})) =cat(1,CELLS.Salsa6fPix{idxC}{:});
Ch3 = Ch3(y_lim,x_lim,:);
him(3) = imshow(Ch3(:,:,idxT),[min(Ch3(:)) max(Ch3(:))],...
    'Parent',ax(9));

hold(ax(7),'on')
hold(ax(8),'on')
hold(ax(9),'on')
h_bnd(1) = line(nan,nan,'Color',[1 1 1],'Parent',ax(7));
h_bnd(2) = copyobj(h_bnd(1),ax(8));
h_bnd(3) = copyobj(h_bnd(1),ax(9));
hold(ax(7),'off')
hold(ax(8),'off')
hold(ax(9),'off')


ax(7).Colormap = [zeros(512,1) [1:512]' zeros(512,1)]/512;
ax(8).Colormap = [[1:512]' zeros(512,1) zeros(512,1)]/512;
ax(9).Colormap = hot(512);
ax(7).CLim = [min(Ch1(:)) max(Ch1(:))];
ax(8).CLim = [min(Ch2(:)) max(Ch2(:))];
ax(9).CLim = [min(Ch3(:)) max(Ch3(:))];

xy = CELLS.Centroid_xy_pix{idxC};
x = nan(metadata.DimT,1);
y = nan(metadata.DimT,1);
x(CELLS.TimeFrame{idxC}) = xy(:,1)*metadata.ResX;
y(CELLS.TimeFrame{idxC}) = xy(:,2)*metadata.ResY;
x = x-mean(x,'omitnan');
y = y-mean(y,'omitnan');
x = [x; nan];
y = [y; nan];
hp = patch(x,y,seconds([t ;nan]),'EdgeColor','interp','LineWidth',1,...
    'Parent',ax(10));       
ax(10).XLabel.String = 'x (\mum)';
ax(10).YLabel.String = 'y (\mum)';
ax(10).Box = 'on';
ax(10).YDir = 'reverse';
h_rp(6) = line(x(idxT),y(idxT),...
    'Marker','o','Color',[0 0 0],'Parent',ax(10));


hcb_r = colorbar(ax(7));
hcb_g = colorbar(ax(8));
hcb_gr = colorbar(ax(9));
hcb_d = colorbar(ax(10));
hcb_r.Title.String  = 'Pixel intensity';
hcb_g.Title.String  = 'Pixel intensity';
hcb_gr.Title.String = 'Pixel intensity';
hcb_d.Title.String  = 'time (s)';

ax(7).Title.String  = 'GCaMP6f';
ax(8).Title.String  = 'tdTomato';
ax(9).Title.String  = 'Salsa6f';
ax(10).Title.String = 'Trajectory';
ax(10).DataAspectRatioMode ='manual';
ax(10).DataAspectRatio = [1 1 1];
ax(10).PlotBoxAspectRatio = [1 1 1];
% him(4) = imshow(salsa_im(y_lim,x_lim,idxT),[min(salsa_im(:)) max(salsa_im(:))],...
%     'Parent',ax(10));
for idxT =1:metadata.DimT
    him(1).CData = Ch1(:,:,idxT);
    him(2).CData = Ch2(:,:,idxT);  
    him(3).CData = Ch3(:,:,idxT);  
    isThere = CELLS.TimeFrame{idxC} == idxT;
    if any(isThere)
        xy = CELLS.XYBoundary{idxC}{isThere};
        h_bnd(1).XData = xy(:,1)-x_lim(1)+1;
        h_bnd(1).YData = xy(:,2)-y_lim(1)+1;
        h_bnd(2).XData = h_bnd(1).XData;
        h_bnd(2).YData = h_bnd(1).YData;
        h_bnd(3).XData = h_bnd(1).XData;
        h_bnd(3).YData = h_bnd(1).YData;
    else
        h_bnd(1).XData = nan;
        h_bnd(1).YData = nan;
        h_bnd(2).XData = nan;
        h_bnd(2).YData = nan;
        h_bnd(3).XData = nan;
        h_bnd(3).YData = nan;
    end
    h_rp(1).XData = hl_salsa.XData(idxT);
    h_rp(1).YData = hl_salsa.YData(idxT);
    h_rp(2).XData = hl_area.XData(idxT);
    h_rp(2).YData = hl_area.YData(idxT);
    h_rp(3).XData = hl_circ.XData(idxT);
    h_rp(3).YData = hl_circ.YData(idxT);
    h_rp(4).XData = hl_shrt.XData(idxT);
    h_rp(4).YData = hl_shrt.YData(idxT);
    h_rp(5).XData = hl_speed.XData(idxT);
    h_rp(5).YData = hl_speed.YData(idxT);
    h_rp(6).XData = x(idxT);
    h_rp(6).YData = y(idxT);
   
drawnow

frame = getframe(hfig);
writeVideo(vid_obj,frame);
end
close(vid_obj)
end
%
%
% idxT = 1;
% him = imshow(stack.Ch2(:,:,idxT),[min(stack.Ch2(:)) 0.3*max(stack.Ch2(:))]);
% h_bnd = gobjects(height(CELLS),1);
% cmap = hsv(height(CELLS));
% cmap = cmap(randperm(height(CELLS)),:);
% for idxC = 1:height(CELLS)
%     h_bnd(idxC) = line(nan,nan,'Color',cmap(idxC,:));    
% end
% 
% for idxT = 1:metadata.DimT
%     him.CData = stack.Ch2(:,:,idxT);
%     for idxC = 1:height(CELLS)
%         isThere = CELLS.TimeFrame{idxC} == idxT;
%         if any(isThere)
%             0
%             xy = CELLS.XYBoundary{idxC}{isThere};
%             h_bnd(idxC).XData = xy(:,1);
%             h_bnd(idxC).YData = xy(:,2);
%         else
%             1
%             h_bnd(idxC).XData = nan;
%             h_bnd(idxC).YData = nan;
%         end
%     end
%     drawnow
% end
% close all
% hold on
% plot(CELLS.Solidity{idxC})
% plot(CELLS.Circularity{idxC})
% plot(CELLS.Eccentricity{idxC})
% plot(CELLS.ShapeRatio{idxC})

end

