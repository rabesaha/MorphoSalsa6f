function metadata = getOMEmetadata_stk(omeMeta,metadata,fullfilepath)
%GETOMEMETADATA_STK Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 23-Mar-2023 NL

% a= omeMeta.dumpXML
% fileID = fopen('myfile.txt','w');
% nbytes = fprintf(fileID,'%s',a)
% fclose(fileID);
% 
[FilePath,FileName,FileExtension] = fileparts(fullfilepath);

NBImage = omeMeta.getImageCount();

for idxImage = 1:NBImage
    metadata(idxImage).FilePath = FilePath;
    metadata(idxImage).FileName = FileName;
    metadata(idxImage).FileExtension = FileExtension;    
    metadata(idxImage).AcquisitionDate = datetime(omeMeta.getImageAcquisitionDate(0).getValue);
    
    metadata(idxImage).Name     = char(omeMeta.getImageName(idxImage-1));
    if isempty(metadata(idxImage).Name)
        metadata(idxImage).Name = metadata(idxImage).FileName;
    end
    
    metadata(idxImage).DimX = omeMeta.getPixelsSizeX(idxImage-1).getValue();
    metadata(idxImage).DimY = omeMeta.getPixelsSizeY(idxImage-1).getValue();
    metadata(idxImage).DimZ = omeMeta.getPixelsSizeZ(idxImage-1).getValue();
    metadata(idxImage).DimC = omeMeta.getPixelsSizeC(idxImage-1).getValue();
    metadata(idxImage).DimT = omeMeta.getPixelsSizeT(idxImage-1).getValue();

    metadata(idxImage).PixelType = char(omeMeta.getPixelsType(idxImage-1));
    metadata(idxImage).DimensionOrder = char(omeMeta.getPixelsDimensionOrder(idxImage-1));

%     if metadata(idxImage).DimZ>1
%         metadata(idxImage).UnitZ = char(omeMeta.getPixelsPhysicalSizeZ(idxImage-1).unit.getSymbol);
%         metadata(idxImage).ResZ  = omeMeta.getPixelsPhysicalSizeZ(idxImage-1).value;
%     end
    metadata(idxImage).TimeUnit    = char(omeMeta.getPlaneDeltaT(idxImage-1,0).unit.getSymbol);

    NBFrame = metadata(idxImage).DimC*metadata(idxImage).DimZ*metadata(idxImage).DimT;
    t = zeros(NBFrame,1);
    z = zeros(NBFrame,1);
    zidx = zeros(NBFrame,1);
    tidx = zeros(NBFrame,1);    
    for idxF = 1:NBFrame
        t(idxF) = omeMeta.getPlaneDeltaT(idxImage-1,idxF-1).value;
%         z(idxF) = omeMeta.getPlanePositionZ(idxImage-1,idxF-1).value;
        zidx(idxF) = omeMeta.getPlaneTheZ(idxImage-1,idxF-1).getValue()+1;
        tidx(idxF) = omeMeta.getPlaneTheT(idxImage-1,idxF-1).getValue()+1;
    end
    t = t-t(1);
    z = z-z(1);
    switch metadata(idxImage).TimeUnit
        case 's'
            metadata(idxImage).time = seconds(t);
    end
    metadata(idxImage).z = z;
    if metadata(idxImage).DimT>1      
        metadata(idxImage).ResT = mean(diff(metadata(idxImage).time));
    end


    metadata(idxImage).ResX = double(omeMeta.getPixelsPhysicalSizeX(idxImage-1).value);
    metadata(idxImage).ResY = double(omeMeta.getPixelsPhysicalSizeY(idxImage-1).value);
    metadata(idxImage).UnitX = char(omeMeta.getPixelsPhysicalSizeX(idxImage-1).unit.getSymbol);
    metadata(idxImage).UnitY = char(omeMeta.getPixelsPhysicalSizeY(idxImage-1).unit.getSymbol);
  

    for idxC = 1:metadata(idxImage).DimC
        metadata(idxImage).ChannelName{idxC} = char(omeMeta.getChannelName(idxImage-1,idxC-1));
    end



end

end

