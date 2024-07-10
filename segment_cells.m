function Lbl = segment_cells(metadata,stack,DefaultOptions)
%SEGMENT_CELLS Summary of this function goes here
%   Detailed explanation goes here
%
%   Nicolas Liaudet
%   Bioimaging Core Facility - UNIGE
%   https://www.unige.ch/medecine/bioimaging/en/bioimaging-core-facility/
% 
%   CC BY-NC 4.0
%
%   v1.0 24-Mar-2023 NL

idxC = find(ismember(metadata.ChannelName,DefaultOptions.ChannelSeg));
im = stack.(['Ch' num2str(idxC)]);
parg = pyargs('model_type', 'cyto',...  
    'gpu',true,...   
    'net_avg',true);
%  'diam_mean',30,...   

pargE = pyargs('channels',[0,0],...
     'do_3D',false,...
      'stitch_threshold',0.10,...
      'z_axis',2,...
      'tile',false);

im = squeeze(im);
mdl = py.cellpose.models.Cellpose(parg);
pyim = py.numpy.array(im,dtype=py.numpy.uint16);
OUT = mdl.eval(pyim,pargE);
Lbl = double(OUT{1});
Lbl = permute(Lbl,[2 3 1]);

end

