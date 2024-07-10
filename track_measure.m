function [CELLS,LBL_FijiROI] = track_measure(metadata,stack,Lbl,FijiROI,DefaultOptions)
%TRACK_MEASURE Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 24-Mar-2023 NL
%
%   v1.1 17-Apr-2023 NL
%   TrackID are in the same order than the FIJI Roi
%
%   v1.2 18-May-2023 NL
%   loop on the actual lbl idx given by FiJi when keeping a single cell per
%   FiJi ROI
%   0based index from Fiji to 1based for Matlab  (ROI pixel coordinates)
%
%   v1.3 25-May-2023 NL
%   Fiji ROI boundaries limited to the image size...
%
%   v1.4 05-Sep-2023 NL
%   CELLS Track ID are not reordered thus now correct matchin with the
%   Fiji'so nes...

IndivCells = cell(metadata.DimT,1);
parfor idxT = 1:metadata.DimT
    tmp = regionprops('table',Lbl(:,:,idxT),'Centroid','Area','Circularity','Solidity','Eccentricity','MajorAxisLength','MinorAxisLength','PixelIdxList');
    tmp.ShapeRatio = tmp.MinorAxisLength./tmp.MajorAxisLength;
    tmp(tmp.Area ==0,:)=[];%remove the "empty label" as lbl is not strictly monotonic, regiopros counts the missing label values
    tmp.TimeFrame = repmat(idxT,[height(tmp) 1]);
    if isempty(tmp)
        continue
    end
    PixelIdxList = cat(1,tmp.PixelIdxList{:});
    PixelIdxList = PixelIdxList+(idxT-1)*metadata.DimX*metadata.DimY;
    PixelIdxList = mat2cell(PixelIdxList, cellfun(@numel,tmp.PixelIdxList),1);
    tmp.PixelIdxList = PixelIdxList;
    IndivCells{idxT} = tmp;
end
IndivCells = cat(1,IndivCells{:});
IndivCells.Properties.VariableNames{ismember(IndivCells.Properties.VariableNames,'Centroid')} = 'Centroid_xy_pix';

% Trackparam.distPixMax = 20;
% Trackparam.mem   = 3;
% Trackparam.good  = 1;
% Trackparam.dim   = 2;
% Trackparam.quiet = 0;
DefaultOptions.Trackparam.mem = 60;
trk = track([IndivCells.Centroid_xy_pix IndivCells.TimeFrame],DefaultOptions.Trackparam.distPixMax,DefaultOptions.Trackparam);

%% Give the track_id to the right object based on their centroids
% Condense the y and x coordinates into a single number
offset = 10^8;
n_trk = trk(:,1)*offset+trk(:,2);
n_Centroid = IndivCells.Centroid_xy_pix(:,1)*offset+IndivCells.Centroid_xy_pix(:,2);
[isFound,idx] = ismember(n_Centroid,n_trk);
if any(~isFound)
    disp('Something is wrong!')
end

IndivCells.TrackID = trk(idx,end);
VariableNames = IndivCells.Properties.VariableNames;
CELLS = groupsummary(IndivCells,'TrackID',@(x) {x});
CELLS.Properties.VariableNames(2:end) = ['Duration_Frame' VariableNames(1:end-1)]; 


%% Filter cells touching the FOV borders
if DefaultOptions.RemovedIfTouchBorder
    FrameTop   = [[1:metadata.DimX]' ones(metadata.DimX,1)];
    FrameBot   = [[1:metadata.DimX]' metadata.DimX*ones(metadata.DimX,1)];
    FrameLeft  = [ones(metadata.DimY,1) [1:metadata.DimY]'];
    FrameRight = [metadata.DimY*ones(metadata.DimY,1) [1:metadata.DimY]'];
    FramePidx  = [FrameTop;FrameBot;FrameLeft;FrameRight];
    FramePidx  = sub2ind([metadata.DimY metadata.DimX],FramePidx(:,2),FramePidx(:,1));
    idxDiscardCell = false(height(CELLS),1);
    for idxC = 1:height(CELLS)
        plist = unique(cat(1,CELLS.PixelIdxList{idxC}{:}));
        if any(ismember(plist,FramePidx))
            idxDiscardCell(idxC) = true; 
        end
    end
    discardedCELLS = CELLS(idxDiscardCell,:);
    CELLS = CELLS(~idxDiscardCell,:);
end

%% Keep CELLS inside Fiji ROI

LBL_FijiROI = zeros(metadata.DimY,metadata.DimX,metadata.DimT,'uint8');
for idxR = 1:length(FijiROI)
    ymin = max(FijiROI{idxR}.vnRectBounds(2)+1,1);
    xmin = max(FijiROI{idxR}.vnRectBounds(1)+1,1);
    ymax = min(FijiROI{idxR}.vnRectBounds(4)+1,metadata.DimY);
    xmax = min(FijiROI{idxR}.vnRectBounds(3)+1,metadata.DimX);
    LBL_FijiROI(xmin:xmax,ymin:ymax,:) = idxR;
end
% isCellKept = false(height(CELLS),1);
InLabel = nan(height(CELLS),1);
for idxC = 1:height(CELLS)
    idxin = cat(1,CELLS.PixelIdxList{idxC}{:});
    val = LBL_FijiROI(idxin);   
    p_in = sum(val~=0)/(length(idxin));
    if p_in>0.5
%         isCellKept(idxC) = true;        
        InLabel(idxC) = mode(val(val~=0));
    end
end
LBL_FijiROI = LBL_FijiROI(:,:,1);



% a = mean(stack.Ch2,4);
% imshow(labeloverlay(imadjust(uint16(a)),max(BW_FijiROI,[],3)),[])
% cmap = hsv(height(CELLS));
% for idxC = 1:height(CELLS)
%     if isCellKept(idxC)
%     line(CELLS.Centroid_xy_pix{idxC}(:,1),CELLS.Centroid_xy_pix{idxC}(:,2),...
%         'Color',cmap(idxC,:))
%     else
%         line(CELLS.Centroid_xy_pix{idxC}(:,1),CELLS.Centroid_xy_pix{idxC}(:,2),...
%         'Color',[1 1 1],'LineStyle',':')
%     end
% end


CELLS = CELLS(~isnan(InLabel),:);

%% Check for discountinuous tracking and merge if possible
% for idxR = 1:max(InLabel)
%     idxR
%     idxC = find(InLabel == idxR);
%     CELLS.TimeFrame(idxC)
% end

%% Keep a single tracked cell per Fiji Roi
InLabel = InLabel(~isnan(InLabel));
LblIdx = unique(InLabel);
idxCellKept = false(height(CELLS),1);

for idxR = 1:length(LblIdx)
    idxC = find(InLabel == LblIdx(idxR));
    [~,idx] = max(cellfun(@length,CELLS.TimeFrame(idxC)));
    % idxCellKept(LblIdx(idxR)) = idxC(idx);
    idxCellKept(idxC(idx)) = true;
end
CELLS = CELLS(idxCellKept,:);

CELLS.TrackID = InLabel(idxCellKept);



%% Get the boundaries 
XYBoundary = cell(height(CELLS),1);
for idxC = 1:height(CELLS)
     BW = false(metadata.DimY,metadata.DimX,metadata.DimT);
     BW(cat(1,CELLS.PixelIdxList{idxC}{:})) = true;
     BW = BW(:,:,CELLS.TimeFrame{idxC});
     ThisXYBoundary = cell(size(BW,3),1); 
     parfor idxT = 1:size(BW,3)
         B = bwboundaries(BW(:,:,idxT),'noholes');         
         ThisXYBoundary{idxT} = [B{1}(:,2) B{1}(:,1)];
    end   
    XYBoundary{idxC}  = ThisXYBoundary;   
end
CELLS.XYBoundary = XYBoundary;

%% Get the intensities
GCaMP6fPix  = cell(height(CELLS),1);
tdTomatoPix = cell(height(CELLS),1);
for idxC = 1:height(CELLS)    
    plist = cat(1,CELLS.PixelIdxList{idxC}{:});
    sz = cellfun(@length,CELLS.PixelIdxList{idxC});    
    for idxCh = 1:metadata.DimC
        im = stack.(['Ch' num2str(idxCh)]);
        val = im(plist);
        val = mat2cell(val,sz);
        switch metadata.ChannelName{idxCh}
            case 'GFP'
                GCaMP6fPix{idxC} = val;
            case 'RFP'
                tdTomatoPix{idxC} = val;
        end

    end
end
CELLS.GCaMP6fPix  = GCaMP6fPix;
CELLS.tdTomatoPix = tdTomatoPix;
Salsa6fPix = cell(height(CELLS),1);
Salsa6fSig = cell(height(CELLS),1);
for idxC = 1:height(CELLS)
    g = CELLS.GCaMP6fPix{idxC};
    r = CELLS.tdTomatoPix{idxC};
    Salsa6fPix{idxC} = cellfun(@(x,y) single(x)./single(y) ,g,r,'UniformOutput',false);    
    Salsa6fSig{idxC} = cellfun(@mean,Salsa6fPix{idxC});
end
CELLS.Salsa6fPix  = Salsa6fPix;
CELLS.Salsa6fSig  = Salsa6fSig;

%% speed, MSD, etc
displacement = cell(height(CELLS),1);
speed = cell(height(CELLS),1);
MSD = cell(height(CELLS),1);
msd_fit = cell(height(CELLS),1);
msd_gof = cell(height(CELLS),1);
lag = cell(height(CELLS),1);
 t = seconds(metadata.time);
 Tclip = 0.25;
for idxC = 1:height(CELLS)
    x = nan(metadata.DimT,1);
    y = nan(metadata.DimT,1);
    x(CELLS.TimeFrame{idxC}) = CELLS.Centroid_xy_pix{idxC}(:,2);
    y(CELLS.TimeFrame{idxC}) = CELLS.Centroid_xy_pix{idxC}(:,2);
    x = x*metadata.ResX;
    y = y*metadata.ResX;

    displacement{idxC} = [NaN; (sum(diff([x y],1,1).^2,2)).^0.5];
    speed{idxC} = displacement{idxC}/seconds(metadata.ResT);

   

    Dt  = zeros(1,length(t)-1);
    msd = zeros(1,length(t)-1);
    ssd = zeros(1,length(t)-1);
    semsd = zeros(1,length(t)-1);
    for dt = 1:length(t)-1
        Dt(dt) = t(1+dt)-t(1);
        dx  = x(1:end-dt)-x(1+dt:end);
        dy  = y(1:end-dt)-y(1+dt:end);
        dr2 = dx.^2+dy.^2;
        msd(dt) = mean(dr2,'omitnan');
        ssd(dt) = std(dr2,'omitnan');

        n = length(t)-dt;
        semsd(dt) =  ssd(dt)/n^0.5;

    end
    MSD{idxC} = msd';
    lag{idxC} = Dt';


    lag_f = Dt(2:round(length(Dt)*Tclip));
    MSD_f = msd(2:round(length(Dt)*Tclip));

    idxKeep = ~isnan(MSD_f);
    lag_f = lag_f(idxKeep);
    MSD_f = MSD_f(idxKeep);
    

    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';%'final' 'iter'
    opts.TolFun = 10^-10;
    opts.TolX   = 10^-10;
    opts.MaxFunEvals = 10^4;
    opts.MaxIter = 10^4;

    diff_model = fittype( '4*D*t.^alpha',...
        'independent', 't', 'dependent', 'msd',...
        'coefficients',{'D','alpha'});
    opts.StartPoint =[ 3e-3   1   ];
%     opts.Lower =     [    0   -50   MSD_f(1)];
%     opts.Upper	=    [   10  50        1e3];

    [msd_fit{idxC}, msd_gof{idxC}] = fit(lag_f', MSD_f', diff_model,opts);


end
CELLS.displacement = displacement;
CELLS.speed = speed;
CELLS.MSD = MSD;
CELLS.lag = lag;
CELLS.MSD_fit = msd_fit;
CELLS.MSD_gof = msd_gof;
%% reorder ID

% CELLS.TrackID = [1:height(CELLS)]';


% %%
% close all
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
