function [fg, keep]=  bsc_tractByEndpointROIs(fg, rois)
% [fg, keep]=  bsc_tractByEndpointROIs(fg, rois)
%
% Purpose:  this function segments streamlines from a whole brain
% fg/tractogram such that a streamline endpoint is in both sets of the
% input rois.  Prevents odd situations that were arising with other
% segmentation methods (i.e. inclusion of within-roi u fibers).
%
% Inputs:
%
% fg:  the tractogram that is to be segmented from.  Typically a whole
% brain fg/tractogram
%
% rois: a pair of rois corresponding to the desired endpoint terminations.
% Input as [{ROI1} {ROI2}]
%
% Outputs:
%
% fg: the output fg structure featuring streamlines meeting the input
% criteria
%
% keep:  a boolean vector locating the indicies in the input fg structure
% of streamlines meeting the criteria
%
%  (C) Daniel Bullock 2018 Bloomington
%% begin code
% hardcode min distance for intersection
minDist = 0.87;

%cut streamlines to just the first and last node
for istreamlines=1:length(fg.fibers)
    endpoint1(:,istreamlines)=fg.fibers{istreamlines}(:,1) ;
    endpoint2(:,istreamlines)=fg.fibers{istreamlines}(:,end);
end

% compute endpoint distnaces to each ROI
[~, distr1e1]=nearpoints(endpoint1, rois{1}.coords');
[~, distr1e2]=nearpoints(endpoint2, rois{1}.coords');

[~, distr2e1]=nearpoints(endpoint1, rois{2}.coords');
[~, distr2e2]=nearpoints(endpoint2, rois{2}.coords');

% only keep those streamlines what have one endpoint meeting each criteria.
keep=or(and(distr1e1<minDist,distr2e2<minDist),and(distr2e1<minDist,distr1e2<minDist));

%transfer relevant streamlines in to output fg
fg.fibers=fg.fibers(keep);
end
