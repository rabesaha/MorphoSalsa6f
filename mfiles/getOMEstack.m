function stack = getOMEstack(metadata,fullfilepath)
%GETOMESTACK Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%
%   v1.0 23-Mar-2023 NL


NBImage = length(metadata);

for idxImage = 1:NBImage
    for idxC=1:metadata(idxImage).DimC
        ch = zeros(metadata(idxImage).DimY,...
            metadata(idxImage).DimX,...
            metadata(idxImage).DimT,...
            metadata(idxImage).DimZ,...
            metadata(idxImage).PixelType);

        for idxZ=1:metadata(idxImage).DimZ

            parfor idxT=1:metadata(idxImage).DimT
                 bfInitLogging('INFO');
                 r2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
                 r2.setId(fullfilepath);
                 r2.setSeries(idxImage-1)
                iplane = loci.formats.FormatTools.getIndex(r2,...
                    idxZ-1,idxC-1,idxT-1)+1;
                ch(:,:,idxT,idxZ) = bfGetPlane(r2, iplane);
                r2.close()
            end
        end
        stack.(['Ch' num2str(idxC)]) = ch;
    end

end

end

