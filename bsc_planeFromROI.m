function [roiOUT] = bsc_planeFromROI(roiIN, location,fsDir)
% DESCRIPTION:
% This function creates a plane using a specified component of an input ROI
% as 
%
% INPUTS:
% -fsDir: path to THIS SUBJECT'S freesurfer directory
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

labelNifti= wma_getAsegFile(fsDir , '2009');


if isnumeric(roiIN)
    
    roiIN=bsc_roiFromFSnums(fsDir,[roiIN ],0);
    
elseif isstruct(roiIN)
    
    %no action
else
    error('Unrecognized input for "roiIN" in bsc_planeFromROI')
end

LRflag=mean(roiIN.coords(:,1))<0;

switch lower(location)
    case 'top'
        roiCoord=max(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'superior'
        roiCoord=max(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'bottom'
        roiCoord=min(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'inferior'
        roiCoord=min(roiIN.coords(:,3));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'z');
    case 'anterior'
        roiCoord=max(roiIN.coords(:,2));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'y');
    case 'posterior'
        roiCoord=min(roiIN.coords(:,2));
        roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'y');
    case 'medial'
        if LRflag
            roiCoord=max(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        else
            roiCoord=min(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        end
    case 'lateral'
        if LRflag
            roiCoord=min(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        else
            roiCoord=max(roiIN.coords(:,1));
            roiOUT=bsc_makePlanarROI(labelNifti,roiCoord, 'x');
        end
        
end

end
            