function [roiOUT] = bsc_planeFromROI_v2(roiIN, location,atlas)
% DESCRIPTION:
% This function creates a plane using a specified component of an input ROI
% as 
%
% INPUTS:
% -atlas:  path to the atlas nii.gz that you would like to use.
%
% -roiIN: the roi that is to serve as reference of the planar ROI.
%
% -location: a string input indicating the extreme-most part of the roi
% that is to serve as the basis for forming the roi.  i.e. "top",
% "superior", "bottom", "inferior", "anterior", "posterior", "medial",
% "lateral"
%
% OUTPUTS:
% -roiOUT: the roi structure of the planar ROI
%
%  (C) Daniel Bullock 2017 Bloomington
%% Begin code
if or(ischar(atlas),isstr(atlas))
  atlas=niftiRead(atlas);
end

if isnumeric(roiIN)
    
    roiIN=bsc_roiFromAtlasNums(atlas,[roiIN ],1);
    
elseif isstruct(roiIN)
    
    %no action
else
    error('Unrecognized input for "roiIN" in bsc_planeFromROI')
end

LRflag=mean(roiIN.coords(:,1))<0;

switch lower(location)
    case 'top'
        roiCoord=max(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'z');
    case 'superior'
        roiCoord=max(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'z');
    case 'bottom'
        roiCoord=min(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'z');
    case 'inferior'
        roiCoord=min(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'z');
    case 'anterior'
        roiCoord=max(roiIN.coords(:,2));
        roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'y');
    case 'posterior'
        roiCoord=min(roiIN.coords(:,2));
        roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'y');
    case 'medial'
        if LRflag
            roiCoord=max(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'x');
        else
            roiCoord=min(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'x');
        end
    case 'lateral'
        if LRflag
            roiCoord=min(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'x');
        else
            roiCoord=max(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(atlas,roiCoord, 'x');
        end
        
end

end
            