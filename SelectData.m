function FullFilePath = SelectData(DefaultOptions)
%SELECTDATA Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 23-Mar-2023 NL
%
%   v1.1 17-Apr-2023 NL
%   files in a "rubbish" folder are now discarded

MainDir = uigetdir(DefaultOptions.LastFolderPath);
if isnumeric(MainDir)
    return
end
fname = dir(fullfile(MainDir,'**','*.stk'));
FileNames = {fname.name}';
FilePaths = {fname.folder}';

if rem(length(FileNames),2) ~=0
    opts = struct('WindowStyle','modal',...
        'Interpreter','none');
    f = errordlg({'No file selected!'},...
        'Analysis stop', opts);
    FullFilePath = [];
    return
end
idxKeep = contains(FileNames,'RFP')|contains(FileNames,'mCherry');
FileNames = FileNames(idxKeep);
FilePaths = FilePaths(idxKeep);
% FileNames = FileNames(1:2:end);
% FilePaths = FilePaths(1:2:end);



DefaultOptions.LastFolderPath = MainDir;
save(['mfiles' filesep 'DefaultOptions.mat'],'DefaultOptions')

FullFilePath = cell(length(FilePaths),1);
for idxF = 1:length(FilePaths)
    FullFilePath{idxF} = fullfile(FilePaths{idxF},FileNames{idxF});
end

FullFilePath(contains(FullFilePath,'rubbish','IgnoreCase',true)) = [];
if isempty(FullFilePath)
    opts = struct('WindowStyle','modal',...
        'Interpreter','none');
    f = errordlg({'No file selected!'},...
        'Analysis stop', opts);
    FullFilePath = [];
    return
end

end

