function [mergedROI] =bsc_roiFromAtlasNums(atlas,ROInums, smoothKernel)
%
%  INPUTS:
%  -fsDir:  path to THIS SUBJECT'S freesurfer rectory.
%
%  -fsROInums:  a list of region ID numbers from freesurfer that you would
%  like merged into a single ROI.  i.e. [ 3, 4 , ... 41,42]
%
%  OUTPUTS:
%  -mergedROI: a merged ROI which has coordinates corresponding to each of
%  the fs regions associated with a number entered in the fsROInums.  For
%  example [3,42] would create a grey matter mask roi, while [2, 41] would
%  create a white matter mask roi
%
%  NOTE:  Makes a call to mri_convert, so requires that FreeSurfer be
%  installed and set up properly.
% % (C) Daniel Bullock 2017 Bloomington, Indiana
%
%% preliminaries



%% set up aparcAsegFile

if or(isstring(atlas),ischar(atlas))
    [fpath,fname,EXT] = fileparts(atlas);
    if strcmp(EXT,'.mgz')
        fprintf('\n .mgz file entered.  Creating .nii.gz using mri_convert')
        %apaprently necessary for matlab?
        spaceChar={' '};
        %I dont know why this is necessary
        quoteString=strcat('mri_convert',spaceChar, atlas ,spaceChar, fullfile(fpath,fname,'.nii.gz'));
        quoteString=quoteString{:};
        [status result] = system(quoteString, '-echo');
        if status~=0
            warning('/n Error generating aseg nifti file.  There may be a problem finding the file. Output: %s ',result)
            
        end
    end
    atlas=niftiRead(atlas);
    fprintf('\n atlas loaded')
else
    %do nothing
end





%% begin loop
fprintf('\n Generating composite roi from region(s) %s', num2str(ROInums))
%get size of atlasNifti.data and make a blank matrix mask for it
atlasDataSize=size(atlas.data);
blankLabelNifti(1:atlasDataSize(1),1:atlasDataSize(2),1:atlasDataSize(3))=false;

ROImask=blankLabelNifti;
roiNameString=[];
for iRois = ROInums
    ROImask=or(ROImask,atlas.data==iRois);
    if iRois==ROInums(end)
        roiNameString=strcat(roiNameString,num2str(iRois));
    else
        roiNameString=strcat(roiNameString,num2str(iRois),'_');
    end
end

%% converts nifti formatting to roi formatting
%smooth if you want to
if ~isempty(smoothKernel) & ~smoothKernel==0
ROImask=~(smooth3(ROImask,'box', smoothKernel))==0;
end

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
