function  [newROI] = bsc_subtractROIs(ROI1, ROI2)
% [roiIntersection] = bsc_intersectROIs(ROI1, ROI2)
%
% Subtract ROI 1 from ROI 2
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
newROI=ROI1;
roiIntersection.coords=[];
roiIntersection.name=strcat(ROI1.name,'_to_',ROI2.name);
if ~isempty(ROI1.coords)
intersectionCoords=intersect(ROI1.coords,ROI2.coords,'rows');
roiIntersection.coords=intersectionCoords;
newROI.coords=setdiff(ROI2.coords,roiIntersection.coords,'rows');
else
    warning('ROI1 is empty')
    newROI=ROI2;
end
end



