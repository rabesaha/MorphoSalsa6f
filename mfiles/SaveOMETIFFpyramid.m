function SaveOMETIFFpyramid(im,Meta,ImageFullPath)
%SAVEOMETIFF Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   v1.0 25-Aug-2021 NL
%
%   v2.0 01-Dec-2022 NL
%   dimension order of im: YXCZT
%
%   v3.0 05-Jan-2023 NL
%   Now work with any dimensions



if exist(ImageFullPath,'file') == 2
    delete(ImageFullPath)
end

%Generate pyramid           
sz = [2 4 8];
I = cell(length(sz)+1,1);
I{1} = im;

for idxP = 1:length(sz)    
    tmp = imresize(im,1/sz(idxP),'Antialiasing',true);%only resize the first 2 dimesions
    I{idxP+1} = tmp;
end


% Create Writer object from output path
writer = loci.formats.out.PyramidOMETiffWriter();



metadata = createMinimalOMEXMLMetadata(I{1},Meta.DimensionOrder);

toInt = @(x) javaObject('ome.xml.model.primitives.PositiveInteger', ...
                        javaObject('java.lang.Integer', x));
for idxP = 2:length(I)
    metadata.setResolutionSizeX(toInt(size(I{idxP},2)), 0, idxP-1);
    metadata.setResolutionSizeY(toInt(size(I{idxP},1)), 0, idxP-1);
    
end
% 
pixelSizeX = ome.units.quantity.Length(java.lang.Double(Meta.ResX),...
    ome.units.UNITS.MICROMETER);
pixelSizeY = ome.units.quantity.Length(java.lang.Double(Meta.ResY),...
    ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeX(pixelSizeX, 0);
metadata.setPixelsPhysicalSizeY(pixelSizeY, 0);

if isfield(Meta,'cmap')
    cmap = Meta.cmap;
else
cmap = [
    31,120,180;...
    51,160,44;...
     227,26,28;...
      255,127,0;...
      106,61,154;...
    178,223,138;...    
    251,154,153;...
    166,206,227;...   
    253,191,111;...   
    202,178,214;...    
    255,255,153;...
    177,89,40];
end



for idxC = 1:Meta.DimC
    metadata.setChannelName(java.lang.String(Meta.ChannelName{idxC}),0,idxC-1);
    jColor = ome.xml.model.primitives.Color(cmap(idxC,1),cmap(idxC,2),cmap(idxC,3),1);    
    metadata.setChannelColor(jColor,0,idxC-1);
end


writer.setWriteSequentially(true);
writer.setMetadataRetrieve(metadata);
writer.setCompression('zlib');
writer.setBigTiff(true);
writer.setId(ImageFullPath);



% Load conversion tools for saving planes
switch class(im)
    case {'int8', 'uint8'}
        getBytes = @(x) x(:);
    case {'uint16','int16'}
        getBytes = @(x) javaMethod('shortsToBytes', 'loci.common.DataTools', x(:), 0);
    case {'uint32','int32'}
        getBytes = @(x) javaMethod('intsToBytes', 'loci.common.DataTools', x(:), 0);
    case {'single'}
        getBytes = @(x) javaMethod('floatsToBytes', 'loci.common.DataTools', x(:), 0);
    case 'double'
        getBytes = @(x) javaMethod('doublesToBytes', 'loci.common.DataTools', x(:), 0);
end


TileDimX = 1024;%1024;
TileDimY = 1024;%1024;
writer.setTileSizeX(TileDimX);
writer.setTileSizeY(TileDimY);

 dimCoord = [Meta.DimC Meta.DimZ Meta.DimT];

for idxP = 1:length(I)

    writer.setResolution(idxP-1);
    isBlockWritting = false;
    if prod(size(I{idxP},[1 2]))>2^31-1%2147483647
        isBlockWritting = true;
    end

    % Save planes to the writer
    nPlanes = metadata.getPixelsSizeZ(0).getValue()*...
        metadata.getPixelsSizeC(0).getValue()*...
        metadata.getPixelsSizeT(0).getValue();

   
    for index = 1:nPlanes
        [i, j, k] = ind2sub(dimCoord, index);
        frame = I{idxP}(:,:,i,j,k)';        
        if isBlockWritting
            bim = blockedImage(frame,'BlockSize',[TileDimY TileDimX]);
           
            for idxY = 1: bim.SizeInBlocks(1)
                for idxX = 1: bim.SizeInBlocks(2)

                    [blockdata,blockinfo] = getBlock(bim, [idxY idxX]);
                    blockinfo.Size = size(blockdata);
                    writer.saveBytes(index-1, getBytes(blockdata),...
                        blockinfo.Start(1)-1,blockinfo.Start(2)-1,...
                        blockinfo.Size(1),blockinfo.Size(2));
                end
            end



        else
            writer.saveBytes(index-1, getBytes(frame));
        end
    end

end


writer.close();
end



