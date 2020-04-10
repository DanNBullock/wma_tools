function [roiOUT]=bsc_shiftROI(roiIN, shiftDir, shiftMag)
% [roiOUT]=bsc_modifyROI(roiIN, shiftDir, shiftMag)
% DESCRIPTION:
% Moves an ROI the specified magnitude in the specified direction.
%
% INPUTS:
%
% -roiIN: the roi that is to be shifted.
%
% -shiftDir:  the direction of the desired ROI shift.  i.e. "up",
% "superior", "down", "inferior", "anterior", "posterior", "medial",
% "lateral"
%
% -shiftMag: an integer indicating the magnitude (in mm) of the desired ROI
% shift
%
% OUTPUTS:
% -roiOUT: the roi structure of the shifted ROI
%
%  (C) Daniel Bullock 2018 Bloomington
%% Begin code
%if roiIN is fsNum get roi, else do nothing.
if isnumeric(roiIN)
    roiIN=bsc_roiFromFSnums(fsDir,[roiIN ],0);
elseif isstruct(roiIN)
    %no action
else
    error('\n Unrecognized input for "roiIN" in bsc_planeFromROI')
end

%determine if this roi is on the left or the right.  Will cause
%problems for cross hemispheric ROIS
LRflag=mean(roiIN.coords(:,1))<0;

%create dummy structure
roiOUT=roiIN;

switch lower(shiftDir)
    case 'up'
        roiOUT.coords=roiIN.coords(:,3)+shiftMag;
    case 'superior'
        roiOUT.coords=roiIN.coords(:,3)+shiftMag;
    case 'down'
        roiOUT.coords=roiIN.coords(:,3)-shiftMag;
    case 'inferior'
        roiOUT.coords=roiIN.coords(:,3)-shiftMag;
    case 'anterior'
        roiOUT.coords=roiIN.coords(:,2)+shiftMag;
    case 'posterior'
        roiOUT.coords=roiIN.coords(:,2)-shiftMag;
    case 'medial'
        if LRflag
            roiOUT.coords=roiIN.coords(:,1)+shiftMag;
        else
            roiOUT.coords=roiIN.coords(:,1)-shiftMag;
        end
    case 'lateral'
        if LRflag
            roiOUT.coords=roiIN.coords(:,1)-shiftMag;
        else
            roiOUT.coords=roiIN.coords(:,1)+shiftMag;
        end
end

roiOUT.name=strcat(roiIn.name,'_shift_',num2str(shiftMag),'_',shiftDir);

end