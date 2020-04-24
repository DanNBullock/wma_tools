function [planarROI]=bsc_makePlanarROI_v3(referenceNifti,mmPlane, dimension)
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
%  4/24/2020  edit:  no longer assumes apcp?  Will yell tho...
%% begin code

%still can't detect it the brain is rotated in subject space.
%detect non-orthogonality of inputAtlas
%for these to be violated things would have to be quite wonky
sumNon1Row1=sum(referenceNifti.qto_xyz(1,abs(referenceNifti.qto_xyz(1,1:3))~=1));
sumNon1Row2=sum(referenceNifti.qto_xyz(2,abs(referenceNifti.qto_xyz(2,1:3))~=1));
sumNon1Row3=sum(referenceNifti.qto_xyz(3,abs(referenceNifti.qto_xyz(3,1:3))~=1));
allRotationSum=sum([sumNon1Row1,sumNon1Row2,sumNon1Row3]);
%arbitrary tolerance
if allRotationSum>.001
    error('orthogonality assumption violated, \n creation of planar rois assumes orthogonal orientation')
end

%taken from modified version of vistasoft's dtiRoiNiftiFromMat.m
bb = mrAnatXformCoords(referenceNifti.qto_xyz, [-(size(referenceNifti.data).*referenceNifti.pixdim)/2; (size(referenceNifti.data).*referenceNifti.pixdim)/2-1]);

%change this to change the spacing of roi coordinates
spacingParameter=.5;
%enumerate coordinates for a flat plane plane at the specified coordinate
switch lower(dimension)
    case 'x' 
        stable= mmPlane;
        yRange=min(bb(:,2)):spacingParameter:max(bb(:,2));
        zRange=min(bb(:,3)):spacingParameter:max(bb(:,3)); 

        [outy,outz]=ndgrid(yRange,zRange);
        outycoords=reshape(outy,numel(outy),1);
        outzcoords=reshape(outz,numel(outz),1);
        outxcoords=stable*ones(size(outzcoords));
        roiCoords=horzcat(outxcoords,outycoords,outzcoords);
        
    case 'y'
        stable= mmPlane;
        xRange=min(bb(:,1)):spacingParameter:max(bb(:,1));
        zRange=min(bb(:,3)):spacingParameter:max(bb(:,3)); 

        [outx,outz]=ndgrid(xRange,zRange);
        outxcoords=reshape(outx,numel(outx),1);
        outzcoords=reshape(outz,numel(outz),1);
        outycoords=stable*ones(size(outzcoords));
        roiCoords=horzcat(outxcoords,outycoords,outzcoords);
    case 'z'
        stable= mmPlane;
        xRange=min(bb(:,1)):spacingParameter:max(bb(:,1));
        yRange=min(bb(:,2)):spacingParameter:max(bb(:,2)); 

        [outx,outy]=ndgrid(xRange,yRange);
        outxcoords=reshape(outx,numel(outx),1);
        outycoords=reshape(outy,numel(outy),1);
        outzcoords=stable*ones(size(outycoords));
        roiCoords=horzcat(outxcoords,outycoords,outzcoords);
end

%set roi Name
roiName=strcat(dimension,'_',num2str(mmPlane),'_plane');

%create the roi Structure
planarROI=dtiNewRoi(roiName,'r',roiCoords);

end 