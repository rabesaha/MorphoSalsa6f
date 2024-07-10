function FijiROI = getFijiROI(FolderPath)
%GETFIJIROI Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 24-Mar-2023 NL


fnames  = dir(fullfile(FolderPath,'*analysis.zip'));
if isempty(fnames)   
    error(['There is no *analysis.zip file in ' FolderPath])
    FijiROI = []
    return
end
FijiROI = ReadImageJROI(fullfile(fnames.folder,fnames.name));

end

