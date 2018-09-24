function  [newROI] = bsc_mergeROIs(ROI1, ROI2)
% [newROI] = bsc_mergeROIs(ROI1, ROI2)
%
% add coords from ROI 1 to ROI 2
%
% Inputs:
% -ROI1:  a vistasoft format ROI 
% -ROI2:  a vistasoft format ROI 
%
% Outputs:
% -newROI:  a vistasoft format ROI, with the coords field
% corresponding to unique corrdinates of the union of the two rois.
%
% (C) Daniel Bullock, 2017, Indiana University
newROI.name=strcat(ROI1.name,'_and_',ROI2.name);

newROI.coords=union(ROI1.coords,ROI2.coords,'rows');
end



