function metadata = correct_metadata(metadata,DefaultOptions)
%CORRECT_METADATA Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 24-Mar-2023 NL

if ~isfield(metadata,'isCorrected')
    metadata.ResZ = DefaultOptions.ResZ;
    metadata.ResX = DefaultOptions.ResX;
    metadata.ResY = DefaultOptions.ResY;
    tmp = metadata.DimZ;
    metadata.DimZ = metadata.DimT;
    metadata.DimT = tmp;
    metadata.time = metadata.time(1:metadata.DimT);
    metadata.ResT = mean(diff(metadata.time));
    metadata.isCorrected = true;

    idx = ismember(metadata.ChannelName,'mCherry');
    if any(idx)       
        metadata.ChannelName(idx) = {'RFP'};
    end
end

end

