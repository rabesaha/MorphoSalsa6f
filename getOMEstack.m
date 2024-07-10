function stack = getOMEstack(metadata,fullfilepath);
%GETOMESTACK Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 23-Mar-2023 NL

% reader.close()

NBImage = length(metadata);

for idxImage = 1:NBImage
%     reader.setSeries(idxImage-1)
    %     cnt = 1;
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
                %                 if etime(clock,last_time) >= 0.25
                %                     progress = cnt/(metadata(idxImage).DimC*metadata(idxImage).DimZ*metadata(idxImage).DimT);
                %                     elap = etime(clock,init_time);
                %                     sec_remain = elap*(1/progress-1);
                %                     r_mes = datestr(sec_remain/86400,'HH:MM:SS');
                %                     waitbar(progress,hwb,{['Reading data ',...
                %                         num2str(cnt) '/' num2str(metadata(idxImage).DimC*metadata(idxImage).DimZ*metadata(idxImage).DimT)];...
                %                         ['Remaining time: ' r_mes]})
                %                     last_time = clock;
                %                     %                end
                %                     if getappdata(hwb,'canceling')
                %                         metadata = [];
                %                         delete(hwb)
                %                         return
                %                     end
                %                     cnt = cnt+1;
                r2.close()
            end
        end
        stack.(['Ch' num2str(idxC)]) = ch;
    end

end

end

