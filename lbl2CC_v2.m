function [CELLS] = lbl2CC_v2(lbl,metadata,stack)
%LBL2CC Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 01-Feb-2023 NL




CELLS = cell(metadata.DimT,1);
parfor idxT = 1:metadata.DimT
    tmp = regionprops('table',lbl(:,:,idxT),'Centroid','Area','Circularity','Solidity','MajorAxisLength','MinorAxisLength','PixelIdxList');
    tmp(tmp.Area ==0,:)=[];%remove the "empty label" as lbl is not strictly monotonic, regiopros counts the missing label values
    tmp.TimeFrame = repmat(idxT,[height(tmp) 1]);

    PixelIdxList = cat(1,tmp.PixelIdxList{:});
    PixelIdxList = PixelIdxList+(idxT-1)*metadata.DimX*metadata.DimY;
    PixelIdxList = mat2cell(PixelIdxList, cellfun(@numel,tmp.PixelIdxList),1);
    tmp.PixelIdxList = PixelIdxList;
    CELLS{idxT} = tmp;
end
CELLS = cat(1,CELLS{:});

param.mem   = 3;
param.good  = 1;
param.dim   = 2;
param.quiet = 0;
trk = track([CELLS.Centroid CELLS.TimeFrame],20,param);

%Give the track_id to the right object based on their centroids
tic
% Condense the y and x coordinates into a single number
vmax = max([trk(:);CELLS.Centroid(:)]);
offset = 10^(length(num2str(vmax))+1);
n_trk = trk(:,1)*offset+trk(:,2);
n_Centroid = CELLS.Centroid(:,1)*offset+CELLS.Centroid(:,2);
[isFound,idx] = ismember(n_Centroid,n_trk);
if any(~isFound)
    disp('Something is wrong!')
end
toc

CELLS.TrackID = trk(idx,end);


VariableNames = CELLS.Properties.VariableNames;
CELLS = groupsummary(CELLS,'TrackID',@(x) {x});
CELLS.Properties.VariableNames(2:end) = ['Duration_Frame' VariableNames(1:end-1)]; 





%% labelled image to Pixel indices of connected components
% 
% CCidx = unique(lbl(:));
% CCidx(CCidx==0) = [];
% NBCC = length(CCidx);
% 
% PixelIdxList = cell(NBCC,1);
% 
% % tic %13
% A = reshape(1:numel(lbl),metadata.DimY,metadata.DimX,metadata.DimT);
% A(lbl==0) = [];
% lbl(lbl==0) = [];
% for idxLbl = 1:length(CCidx)
%     k_idx = lbl==CCidx(idxLbl);
%     pidx = A(k_idx);
%     PixelIdxList{idxLbl} = pidx';
%     A(k_idx)   = [];
%     lbl(k_idx) = [];
% end
% toc


%%
allPixelIdxList = cat(1,CELLS.PixelIdxList{:});
allPixelIdxList = cat(1,allPixelIdxList{:});
DimYXT  = [metadata.DimY metadata.DimX metadata.DimT];

%__________________________________________________________________________
%convert rapidly PixelIdxList into YXT coordinates
cumdim  = [1 cumprod(DimYXT(1:end-1))];
allPixelIdxList = allPixelIdxList(:) - 1;
YXT = zeros(length(allPixelIdxList),numel(cumdim));
for idxD = numel(cumdim):-1:2
    r = rem(allPixelIdxList , cumdim(idxD));
    YXT(:,idxD) = (allPixelIdxList-r)/cumdim(idxD);
    allPixelIdxList  = r;
end
YXT(:,1) = allPixelIdxList;
YXT = YXT + 1;
YXT = mat2cell(YXT,cellfun(@numel,cat(1,CELLS.PixelIdxList{:})),3);
YXT = mat2cell(YXT,cellfun(@numel,CELLS.PixelIdxList),1);
%__________________________________________________________________________
%positions and area
NBCC = height(CELLS);
StartTidx  = zeros(NBCC,1);%start time index
StopTidx   = zeros(NBCC,1);%stop time index
TIME       = cell(NBCC,1);%time vector in physical units
DurationTP = zeros(NBCC,1);%Duration in time point
phyCe_yxt  = cell(NBCC,1);%centroids in physical units
pixCe_yxt  = cell(NBCC,1);%centroids in pixel
pixArea    = cell(NBCC,1);%Area (pixel) vs time
phyArea    = cell(NBCC,1);%Area (um2) vs time
for idxCC = 1:NBCC
   
    yxt = YXT{idxCC};
    yxt = cat(1,yxt{:});
    t = yxt(:,3);
     
    StartTidx(idxCC) = t(1);
    StopTidx(idxCC)  = t(end);
    
    DurationTP(idxCC) = t(end)-t(1)+1;
    TIME{idxCC} = metadata.time(StartTidx(idxCC):StopTidx(idxCC));
    pixCe_yxt{idxCC}  = splitapply(@(y,x) [mean(y) mean(x)],yxt(:,1),yxt(:,2),findgroups(yxt(:,end)));
    phyCe_yxt{idxCC}  = pixCe_yxt{idxCC}.*[metadata.ResY metadata.ResX];

    pixArea{idxCC} = groupcounts(t);
    phyArea{idxCC} = pixArea{idxCC}*prod([metadata.ResY metadata.ResX]);
   
end

%__________________________________________________________________________
%intensities
for idxC = 1:metadata.DimC    
    im = stack.(['Ch' num2str(idxC)]);
    val = cell(NBCC,1);
    for idxCC = 1:NBCC
        val{idxCC} = single(im(cat(1,CELLS.PixelIdxList{idxCC}{:})));
    end
    switch metadata.ChannelName{idxC}
        case 'GFP'
            GCAMP6fPix = val;
        case 'RFP'
            tdTomatoPix = val;
    end
end
CalSigPix = cellfun(@(ch1,ch2) ch1./ch2,GCAMP6fPix,tdTomatoPix,'UniformOutput',false);
CalSig = cell(NBCC,1);
GCAMP6fSig  = cell(NBCC,1);
tdTomatoSig = cell(NBCC,1);
for idxCC = 1:NBCC
    yxt = cat(1,YXT{idxCC}{:});    
    idxG = findgroups(yxt(:,end));
    CalSig{idxCC} = splitapply(@(val) mean(val),CalSigPix{idxCC},idxG);
    GCAMP6fSig{idxCC}  = splitapply(@(val) mean(val),GCAMP6fPix{idxCC},idxG);
    tdTomatoSig{idxCC} = splitapply(@(val) mean(val),tdTomatoPix{idxCC},idxG);
end

%__________________________________________________________________________
%speed, MSD, etc

displacement = cell(NBCC,1);
speed = cell(NBCC,1);
MSD = cell(NBCC,1);
lag = cell(NBCC,1);
for idxCC = 1:NBCC
    x = phyCe_yxt{idxCC}(:,2);
    y = phyCe_yxt{idxCC}(:,1);
    displacement{idxCC} = [NaN; (sum(diff([x y],1,1).^2,2)).^0.5];
    
    speed{idxCC} = displacement{idxCC}/seconds(metadata.ResT);

    t = seconds(TIME{idxCC});

    Dt  = zeros(1,length(t)-1);
    msd = zeros(1,length(t)-1);
    ssd = zeros(1,length(t)-1);
    semsd = zeros(1,length(t)-1);
    for dt = 1:length(t)-1
        Dt(dt) = t(1+dt)-t(1);
        dx  = x(1:end-dt)-x(1+dt:end);
        dy  = y(1:end-dt)-y(1+dt:end);
        dr2 = dx.^2+dy.^2;
        msd(dt) = mean(dr2);
        ssd(dt) = std(dr2);

        n = length(t)-dt;
        semsd(dt) =  ssd(dt)/n^0.5;

    end
    MSD{idxCC} = msd';
    lag{idxCC} = seconds(Dt)';

end



ID = [1:NBCC]';


tmp = table(YXT,TIME,...
    StartTidx,StopTidx,DurationTP,...
    pixCe_yxt,phyCe_yxt,...
    pixArea,phyArea,...
    GCAMP6fPix,tdTomatoPix,CalSig,...
    displacement,speed,MSD,lag);

CELLS = [CELLS tmp];

end

