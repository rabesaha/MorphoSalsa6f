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
% I{1} = reshape(im , [size(im,1) size(im,2) 1 size(im,3) 1]);
for idxP = 1:length(sz)    
    tmp = imresize(im,1/sz(idxP),'Antialiasing',true);%only resize the first 2 dimesions
    I{idxP+1} = tmp;%reshape(tmp, [size(tmp,1) size(tmp,2) 1 size(tmp,3) 1]);    
end



% Create Writer object from output path
writer = loci.formats.out.PyramidOMETiffWriter();



metadata = createMinimalOMEXMLMetadata(I{1},Meta.DimensionOrder);
% 
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
            % [blockdata,a] = getBlock(bim,[1 6]);
            for idxY = 1: bim.SizeInBlocks(1)
                for idxX = 1: bim.SizeInBlocks(2)
%                     [idxY idxX]
                    [blockdata,blockinfo] = getBlock(bim, [idxY idxX]);
                    blockinfo.Size = size(blockdata);
                    writer.saveBytes(index-1, getBytes(blockdata),...
                        blockinfo.Start(1)-1,blockinfo.Start(2)-1,...
                        blockinfo.Size(1),blockinfo.Size(2));
                end
            end


%             for idxY = 1:TileDimY:size(frame,1)
%                 for idxX = 1:TileDimX:size(frame,2)
%                     
%                     tile = frame(idxY:min(idxY+TileDimY-1,size(frame,1)-1),...
%                                  idxX:min(idxX+TileDimX-1,size(frame,2)-1) );
%                     
%                     [idxX idxY ,size(tile,2),size(tile,1), size(frame,2),size(frame,1) ]
% %                     imshow(tile,[])
% %                     pause
%                     
% 
%                     writer.saveBytes(index-1, getBytes(tile),...
%                         idxY-1,idxX-1,size(tile,1),size(tile,2));
%                     
%                 end
%             end
% 



        else
            writer.saveBytes(index-1, getBytes(frame));
        end
    end

% 
%     if prod(size(I{idxP},[1 2]))>=
% 
%         for index = 1 : nPlanes
%             [i, j, k] = ind2sub(zctCoord, index);
%             plane = I{idxP}(:,:,i,j,k)';
%             %         index
%             for idxY = 1:TileDimXY(2):size(tmp,1)
%                 for idxX = 1:TileDimXY(1):size(tmp,2)
% 
%                     tmptile = tmp(idxY:min(idxY+TileDimXY(2)-1,size(tmp,1)),...
%                         idxX:min(idxX+TileDimXY(1)-1,size(tmp,2)));
%                     imshow(tmptile,[])
%                     pause
%                     %                           try
%                     [idxX idxY ,size(tmptile,2),size(tmptile,1),size(tmp,1),size(tmp,2) ]
%                     writer.saveBytes(index-1, getBytes(tmptile),...
%                         idxX-1,idxY-1,size(tmptile,2),size(tmptile,1));
%                     %                           catch
%                     %                               16516
%                     %                           end
% 
%                 end
%             end
% 
%         end
% 
% 
% 
% 
% 
%     else
% 
%         for index = 1 : nPlanes
%             [i, j, k] = ind2sub(zctCoord, index);
%             plane = I{idxP}(:,:,i,j,k)';
% 
%             writer.saveBytes(index-1, getBytes(plane));
%         end
% 
%     end
end


writer.close();
end



%%

%SizeX="7168" SizeY="5632"

% clc
% Nbx = floor(size(frame,2)/TileDimX)
% Nby = floor(size(frame,1)/TileDimY)
% disp("X: " + Nbx+" tiles of " + TileDimX + " - 1 tile of " + (size(frame,2)-Nbx*TileDimX))
% disp("Y: " + Nby+" tiles of " + TileDimY + " - 1 tile of " + (size(frame,1)-Nby*TileDimY))
% 
% TileDimX*([1:Nbx]-1)+1
% TileDimY*([1:Nby]-1)+1
% 
% disp('___________')
% for idxY = 1:TileDimY:size(frame,1)
% 
% 
%     for idxX = 1:TileDimX:size(frame,2)
% 
%         tile = frame(idxY:min(idxY+TileDimY-1,size(frame,1)-1),...
%             idxX:min(idxX+TileDimX-1,size(frame,2)-1) );
% 
%         disp([idxX idxY ,size(tile,2),size(tile,1), size(frame,2),size(frame,1) ])
%         %                     imshow(tile,[])
%         %                     pause
% 
% 
%       
% 
%     end
% end
% 
% %%
% 
% bim = blockedImage(frame,'BlockSize',[TileDimY TileDimX]);
% % [blockdata,a] = getBlock(bim,[1 6]);
% for idxY = 1: bim.SizeInBlocks(1)  
%     for idxX = 1: bim.SizeInBlocks(2)  
%         [idxY idxX]
%         [blockdata,a] = getBlock(bim, [idxY idxX]);
%         a.Size = size(blockdata)
%     end
% end