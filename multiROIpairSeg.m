function [classificationOUT]=multiROIpairSeg(feORwbfg,varargin)
% function [classificationOUT]=multiROIpairSeg(feORwbfg,varargin)
%
%  Purpose:  This function segments tracts from an input wbfg by
%  iteratively segmenting between input ROI pairs.  In theory, input ROIs
%  can be either .mat format or nifti format, either as strings or as
%  objects.
%
%  Inputs
%
%  feORwbfg:  either a whole brain fiber group / tractogram
%  path/object or an fe path / object.
%
%  varargin:  a series of either .mat rois or niftis, either as strings or
%  as objects.  Inputs are interpreted such that the first input is the
%  first item  of the pair pairing while the second input is the second
%  item of the pair.  Thus odds are ROI1 and evens are ROI2
%
%  Outputs:
%
%  classificationOUT:  A classification structure with names generated
%  'tract_1' through 'tract_n' where n is the number of ROI pairings
%  entered.  Future versions can generate names more intelligently.
%  WARNING: IF A STREAMLINE IS SEGMENTED INTO MULTIPLE TRACTS
%  IT WILL ONLY BE CLASSIFIED UNDER THE MOST RECENT (LAST) TRACT.
%
% %  (C) Daniel Bullock, 2018, Indiana University
%% begin code
%load and parse input fg
[wbFG, ~] = bsc_LoadAndParseFiberStructure(feORwbfg);
%initialize fg structure
classificationOUT=[];
for iPairs=1:length(varargin)/2
    %obtain and load ROIs
    roi1= bsc_loadAndParseROI(ROIorNiftivarargin(iPairs*2-1));
    roi2=bsc_loadAndParseROI(ROIorNiftivarargin(iPairs*2));
    % do standard segmentation
    [~, keep]=  bsc_tractByEndpointROIs(wbFG, [roi1,roi2]);
    %use output boolean to add to classification structure
    [classificationOUT]=bsc_concatClassificationCriteria(classificationOUT,strcat('tract_',num2str(iPairs)),keep);
end
