function  [roiIntersection] = bsc_MultiIntersectROIs(Atlas,Inflate,varargin)
%[roiIntersection] = bsc_MultiIntersectROIs(Atlas,Inflate,varargin)
%
% this function finds the overlap of multiple ROIs
%
% Inputs:
% Atlas:  nifti parcellation or path to such
% Inflate: desired inflation kernel
% varargin: list of rois or atlas indicies
%
% Outputs:
% -roiIntersection:  a vistasoft format ROI, with the coords field
% corresponding to the overlap of all rois.
%
% (C) Daniel Bullock, 2017, Indiana University


if isnumeric(varargin{1})
    roiIntersection=bsc_roiFromAtlasNums(Atlas,varargin{1}, Inflate);
else
    roiIntersection=varargin{1};
end

for IROIS=2:length(varargin)
    if  isnumeric(varargin{IROIS})
        curRoi=bsc_roiFromAtlasNums(Atlas,varargin{IROIS}, Inflate);
    else
        curRoi=varargin{IROIS};
    end
    
    [roiIntersection] = bsc_intersectROIs(roiIntersection, curRoi);
end