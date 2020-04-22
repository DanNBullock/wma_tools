function [mergedROI] =bsc_roiFromAtlasNums(atlas,ROInums, smoothKernel)
%  [mergedROI] =bsc_roiFromAtlasNums(atlas,ROInums, smoothKernel)
%
%  INPUTS:
%  -atlas:  path to the atlas nii.gz that you would like to use.
%
%  -ROInums:  a list of region ID numbers from the parcellation that you would
%  like merged into a single ROI.  i.e. [ 3, 4 , ... 41,42]
%
%  -smoothKernel:  a smoothing kernel applied to each coordinate of an roi.
%   Results in an inflation.  NOTE: uses the matlab function smooth3, which
%   places smooths via a box gaussian kernel of edge size [smoothKernel],
%   also, must be an odd integer value, as the inflation occurs evenly in
%   all directions (and thus is odd, due to the inclusion of the origin
%   point).
%
%  OUTPUTS:
%  -mergedROI: a merged ROI which has coordinates corresponding to each of
%  the atlas regions from the passed in label indexes
%
%  NOTE:  Previously atlas input was fsDir input, and had a case for
%  converting .mgz.nii.gz.  Now we will presume .nii.gz (path or object) being passed in
%  this field and throw an error otherwise.
% % (C) Daniel Bullock 2017 Bloomington, Indiana
%

%% get the nifti
if or(isstring(atlas),ischar(atlas))
    %fileparts doesn't work too well here. handles either compressed or
    %uncompressed
     niftiBool= or(strcmp(atlas(end-5:end),'nii.gz'),strcmp(atlas(end-2:end),'nii'));
    %if a nifti is detected load it
    if niftiBool
    atlas=niftiRead(atlas);
    %yell at the user for being inefficeint.  Actually the coder.
    warning('pass a pre-loaded nifti to speed up processing')
    fprintf('\n atlas loaded')
    else %not actually a nifti
        error('Atlas path input (%s) does not point to a nifti',atlas)
    end
else
    %do nothing, its already an object
end

%% generate inermediary nifti mask from requested labels
fprintf('Generating composite roi from region(s) %s\n', num2str(ROInums))
%get size of atlasNifti.data and make a blank matrix mask for it
atlasDataSize=size(atlas.data);
blankLabelNifti(1:atlasDataSize(1),1:atlasDataSize(2),1:atlasDataSize(3))=false;

ROImask=blankLabelNifti;
roiNameString=[];
for iRois = ROInums
    ROImask=or(ROImask,atlas.data==iRois);
    %warn the user if that index isn't found in this nifti
    if ~any(atlas.data==iRois)
        warning('no voxels for %i found in this atlas',iRois)
        fprintf('\n fname: %s',atlas.fname)
        fprintf('\n descrip: %s',atlas.descrip)
    end
        
    if iRois==ROInums(end)
        roiNameString=strcat(roiNameString,num2str(iRois));
    else
        roiNameString=strcat(roiNameString,num2str(iRois),'_');
    end
end


%% converts nifti formatting to roi formatting
%smooth if you want to
initROIsize=length(find(ROImask));
if ~isempty(smoothKernel) & ~smoothKernel==0
    ROImask=~(smooth3(ROImask,'box', smoothKernel))==0;
end
increasePct=((length(find(ROImask))/initROIsize)-1)*100;
fprintf('ROI size increased by %4.2f percent\n',increasePct)
%get index coordinates from nifti mask volume
Roi = find(ROImask);
[x1,y1,z1] = ind2sub(size(atlas.data), Roi);

%create new ROI structure
mergedROI = dtiNewRoi(['', roiNameString], 'r');
%apply the qto_xyz transform to transform the mask coordinates to acpc
%cordinates
mergedROI.coords = mrAnatXformCoords(atlas.qto_xyz, [x1,y1,z1]);
if length(mergedROI.coords)==0
    warning('Empty ROI returned for roi %s', roiNameString)
else
end

end
