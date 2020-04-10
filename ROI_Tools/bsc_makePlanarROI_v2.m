function [planarROI]=bsc_makePlanarROI_v2(fsDir,coord, dimension)
% [planarROI]=bsc_makePlanarROI_v2(fsDir,coord, dimension)
% DESCRIPTION:
% This function creates a plane 
%
% INPUTS:
% -fsDir: path to THIS SUBJECT'S freesurfer directory
%
% -coord: the XYZ ACPC mm coordinate that you would like to generate a planar ROI
% at.  The ROI will be expanded along two of the dimensions associated with
% the plane specified in dimension.  
%
% -dimension: either 'x', 'y', or 'z', to indicate the plane that you would
% like the roi generated along.
%
% OUTPUTS:
% -planarROI: the roi structure of the planar ROI
%
%  (C) Daniel Bullock 2017 Bloomington
%% preliminaries

referenceNifti= wma_getAsegFile(fsDir , '2009');

if ~isstruct(coord)
    sizeLabelNifti=size(referenceNifti.data);

%% make a plane in a volume
% the plane of interest doesnt really matter here, will be replaced later.
switch lower(dimension)
    case 'x'
        blankcheck(1:sizeLabelNifti(2),1:sizeLabelNifti(3))=true;
        mmPlane=coord(1);
    case 'y'
        blankcheck(1:sizeLabelNifti(1),1:sizeLabelNifti(3))=true;
        mmPlane=coord(2);
    case 'z'
        blankcheck(1:sizeLabelNifti(1),1:sizeLabelNifti(2))=true;
        mmPlane=coord(3);
end






blankindex=find(blankcheck);
[x1,y1]=ind2sub(size(blankcheck),blankindex);

%create new xform structure
modXform=zeros(4,4);
modXform(1,1)=referenceNifti.pixdim(1);
modXform(2,2)=referenceNifti.pixdim(2);
modXform(3,3)=referenceNifti.pixdim(3);
%a gamble
modXform(4,4)=referenceNifti.pixdim(3);

modXform(1,4)=-abs(referenceNifti.qoffset_x);
modXform(2,4)=-abs(referenceNifti.qoffset_y);
modXform(3,4)=-abs(referenceNifti.qoffset_z);

% perform transform on coordinates
% and
% set the row  back to the desired value
switch lower(dimension)
    case 'x'
        p1X=zeros(length(x1),1);
        [roiCoords, outMat, imScale] = mrAnatXformCoords(modXform, [p1X, x1, y1]);
        roiCoords(:,1)=mmPlane;
    case 'y'
        p1Y=zeros(length(x1),1);
        [roiCoords, outMat, imScale] = mrAnatXformCoords(modXform, [x1, p1Y, y1]);
        roiCoords(:,2)=mmPlane;
    case 'z'
        p1Z=zeros(length(x1),1);
        [roiCoords, outMat, imScale] = mrAnatXformCoords(modXform, [x1, y1, p1Z]);
        roiCoords(:,3)=mmPlane;
end

%set roi Name
roiName=strcat(dimension,'_',num2str(mmPlane),'_plane');

%create the roi Structure
planarROI=dtiNewRoi(roiName,'r',roiCoords);

end 