function  export_data(metadata,stack,CELLS,LBL_FijiROI)
%EXPORT_DATA Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 27-Mar-2023 NL



% for idxC = 1:height(CELLS)
% end

VarNames = {'Salsa6fSig','Area','Centroid_xy_pix','MajorAxisLength','MinorAxisLength','Eccentricity'...
    'Circularity','Solidity','ShapeRatio','speed'};

for idxV = 1:length(VarNames)
    Meas = CELLS.(VarNames{idxV});

    if strcmp(VarNames{idxV},'Centroid_xy_pix')
        DATA = cell(metadata.DimT+1,2*height(CELLS)+1);
    else
        DATA = cell(metadata.DimT+1,height(CELLS)+1);
    end
    DATA(1,1) = {'time (s)'};
    DATA(2:end,1) = num2cell(seconds(metadata.time));
    for idxC = 1:height(CELLS)
        switch VarNames{idxV}
            case 'Salsa6fSig'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC};
            case 'Area'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC}*metadata.ResX*metadata.ResY;
            case 'Centroid_xy_pix'
                y = nan(metadata.DimT,2);
                y(CELLS.TimeFrame{idxC},:) = Meas{idxC}*metadata.ResX;
            case 'MajorAxisLength'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC}*metadata.ResX;
            case 'MinorAxisLength'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC}*metadata.ResX;
            case 'Eccentricity'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC};
            case 'Circularity'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC};
            case 'Solidity'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC};
            case 'ShapeRatio'
                y = nan(metadata.DimT,1);
                y(CELLS.TimeFrame{idxC}) = Meas{idxC};
            case 'speed'                
                y = Meas{idxC};

        end
        
        if strcmp(VarNames{idxV},'Centroid_xy_pix')
            DATA(1,2*(idxC-1)+1) = {['Cell ' num2str(CELLS.TrackID(idxC),'%02.0f')]};
            DATA(2:end,2*(idxC-1)+[1 2]) = num2cell(y);
        else
            DATA(1,idxC+1) = {['Cell ' num2str(CELLS.TrackID(idxC),'%02.0f')]};
            DATA(2:end,idxC+1) = num2cell(y);
        end

    end
    writecell(DATA,fullfile(metadata.FilePath,['meas_' metadata.FileName '.xlsx']),'Sheet',VarNames{idxV})

   
end

 save(fullfile(metadata.FilePath,[metadata.FileName,'.mat']),...
        'metadata','stack','CELLS','LBL_FijiROI','-v7.3','-nocompression')

end

