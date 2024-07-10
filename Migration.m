%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   Cellpose, https://cellpose.readthedocs.io/en/latest/
%   Bioformat, http://www.openmicroscopy.org/bio-formats/  
%   John C. Crocker, https://physics.emory.edu/faculty/weeks/idl/
%
%   v1.0 23-Mar-2023 NL
%
%   v1.1 17-Apr-2023 NL
%   Modifications in SelectData.m ("rubbish" discard)
%   Modifications in track_measure.m (Fiji ROI order)
%
%   v1.2 18-May-2023 NL
%   loop on the actual lbl idx given by FiJi when keeping a single cell per
%   FiJi ROI
%
%   v1.3 25-May-2023 NL
%   minor modification for mCherry and RFP
%
%   v1.4 08-Aug-2023 NL
%   Takes into account zoom and binning changes
%
%   v1.5 05-Sep-2023 NL
%   TrackID with Fiji ID fixed

%% Initialization
addpath(genpath('mfiles'))
DefaultOptions = Initialization();

%% ------------------- PARAMETERS -------------------


%_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
DefaultOptions.Binning = 2;
DefaultOptions.Zoom = 1.5;
DefaultOptions.RemovedIfTouchBorder = true;
%_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

R_nominal = ((100/240)/4)*1.5; %in um with binning 4 and 1.5 zoom
DefaultOptions.ResX = R_nominal*DefaultOptions.Binning/DefaultOptions.Zoom;
DefaultOptions.ResY = R_nominal*DefaultOptions.Binning/DefaultOptions.Zoom;
DefaultOptions.ResZ = 1; %in um

% DefaultOptions.ChannelSeg = 'mCherry';
DefaultOptions.ChannelSeg = 'RFP';%Do not change
%DefaultOptions.MinTrackDuration = 2/3; % Of the whole acquisition
%duration..... Not used?


%% Select data (stk files) within folders 
FullFilePath = SelectData(DefaultOptions);
if isempty(FullFilePath)
    return
end

%% Main loop over the data...
for idxF = 1:length(FullFilePath)
    disp('_______________________________________________________________')
     disp(['Processing: ',...
            FullFilePath{idxF}])

    try
        tic
        [metadata, stack] = bfread(FullFilePath{idxF});
        toc

        metadata = correct_metadata(metadata,DefaultOptions);

        FijiROI = getFijiROI(fileparts(FullFilePath{idxF}));

        tic
        Lbl = segment_cells(metadata,stack,DefaultOptions);
        toc

        tic
        [CELLS,LBL_FijiROI] = track_measure(metadata,stack,Lbl,FijiROI,DefaultOptions);
        toc

        tic
        display_cells(metadata,stack,CELLS,LBL_FijiROI)
        toc

        tic
        export_data(metadata,stack,CELLS,LBL_FijiROI)
        toc
    catch
        disp(['Something went wrong with: ',...
            FullFilePath{idxF}])

    end


end