function  [roiIntersection] = bsc_intersectROIs(ROI1, ROI2)
% [roiIntersection] = bsc_intersectROIs(ROI1, ROI2)
%
% This function finds the intersection of two ROIs and creates a new ROI
% representing that interection.
%
% Inputs:
% -ROI1:  a vistasoft format ROI 
% -ROI2:  a vistasoft format ROI 
%
% Outputs:
% -roiIntersection:  a vistasoft format ROI, with the coords field
% corresponding to the overlap of the two rois.
%
% (C) Daniel Bullock, 2017, Indiana University

roiIntersection=ROI1;
roiIntersection.coords=[];
roiIntersection.name=strcat(ROI1.name,'_to_',ROI2.name);
intersectionCoords=intersect(ROI1.coords,ROI2.coords,'rows');
roiIntersection.coords=intersectionCoords;
end



