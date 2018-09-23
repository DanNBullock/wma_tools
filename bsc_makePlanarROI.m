function [planarROI]=bsc_makePlanarROI(referenceNifti,mmPlane, dimension)
%
% INPUTS:
% -referenceNifti:  the nifti (or path to nifti) that the ROI will be
% applied to, also functions as the source of affine transform.
%
% -mmPlane: the ACPC mm plane that you would like to generate a planar ROI
% at.  i.e. mmPlane=0 and dimension= x would be a planar roi situated along
% the midsaggital plane.
%
% -dimension: either 'x', 'y', or 'z', to indicate the plane that you would
% like the roi generated along
%
% OUTPUTS:
% -planarROI: the roi structure of the planar ROI
%
%  (C) Daniel Bullock 2017 Bloomington
%% preliminaries
if ischar(referenceNifti)
    referenceNifti = niftiRead(referenceNifti);
else
end



sizeLabelNifti=size(referenceNifti.data);
blankLabelNifti(1:sizeLabelNifti(1),1:sizeLabelNifti(2),1:sizeLabelNifti(3))=false;
%% make a plane in a volume
% the plane of interest doesnt really matter here, will be replaced later.
switch lower(dimension)
    case 'x'
        blankcheck(1:sizeLabelNifti(2),1:sizeLabelNifti(3))=true;
    case 'y'
        blankcheck(1:sizeLabelNifti(1),1:sizeLabelNifti(3))=true;
    case 'z'
        blankcheck(1:sizeLabelNifti(1),1:sizeLabelNifti(2))=true;
       
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