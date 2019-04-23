function [outBool ]=bsc_endpointAtlasCriteria(wbfg,fsDir,roiNums,criteria)
%[outBool ]=bsc_endpointAtlasCriteria(wbfg,fsDir,roiNums,criteria)
%
%  This function determines which streamlines end in the spcecified rois.
%  NOT RECOMMENDED FOR LARGER TRACTS/FGS
%
%  INPUTS
%  wbfg:  an fg structure
%
%  fsDir:  the freesurfer directory corresponding to the same subject (as
%  the source fg).
% 
%  roiNums:  The roi numbers of interest
%
%  criteria:  either 'one', 'both', 'neither', or 'either'
%
%  OUTPUTS:
%
%  outBool:  Boolean corresponding to the streamlines which meat the
%  criteria
%
%  Dan Bullock 2019

[endZone1, endZone2 ]=bsc_investigateTract(wbfg,fsDir);

dualROIVec=[endZone1 endZone2];

appliedCriteria=ismember(dualROIVec,roiNums);

switch lower(criteria)
    case 'one'
        outBool=sum(appliedCriteria,2)==1;
    case 'both'
        outBool=sum(appliedCriteria,2)==2;
    case 'neither'
        outBool=sum(appliedCriteria,2)==0;
    case 'either'
        outBool=~sum(appliedCriteria,2)==0;
end

outBool=outBool;

end
