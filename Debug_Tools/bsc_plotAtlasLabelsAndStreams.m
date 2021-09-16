function bsc_plotAtlasLabelsAndStreams(streamObj,atlas,labels)
% [classificationOut] = bsc_segmentAntPostTracts_v4(wbfg, atlas, categoryClassification)
%
% This function plots atlas label ROIs (as single point cloud) along with input
% streamlines in order to assess the proximety of these objects to one
% another
%
% Inputs:
% -streamObj: a whole brain fiber group structure
% -atlas: path to the desired atlas or the atlas object itself
% -labels: a vector of integers (or a single integer) corresponding to the
% desired Atlas labels to be visualized
%
% Outputs:
% None, produces a plot
%
% (C) Daniel Bullock, 2021, University of Minnesota
%% begin code

%load the fg, if necessary
[fg, ~] = bsc_LoadAndParseFiberStructure(streamObj);

% get the atlas labels as an roi
[mergedROI] =bsc_roiFromAtlasNums(atlas,labels, 1);

% plot the fiber group and the roi
bsc_plotROIandFG(fg,mergedROI,'g')

end

