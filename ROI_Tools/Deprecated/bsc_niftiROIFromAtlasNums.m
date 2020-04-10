function [ni, roiName] =bsc_niftiROIFromAtlasNums(atlas,ROInums, smoothKernel)
%function [mergedROI] =bsc_niftiROIFromAtlasNums(atlas,ROInums, smoothKernel)
%
% Instead of generating a vistasoft compliant roi, generates a a nifti and
% saves it down
%
%  INPUTS:
%  -fsDir: path to an atlas.
%
%  -ROInums:  a list of region ID numbers from the atlas that you would
%  like merged into a single ROI.  i.e. [ 3, 4 , ... 41,42].  Must
%  correspond to the values in the atlas nifti itself. 
%
%  OUTPUTS: -mergedROI: a merged ROI which has coordinates corresponding to
%  each of the fs regions associated with a number entered in the ROInums.
%  For example in the free surfer atlas [3,42] would create a grey matter
%  mask roi, while [2, 41] would create a white matter mask roi.
%
%  NOTE:  Makes a call to mri_convert, so requires that FreeSurfer be
%  installed and set up properly.

%  (C) Daniel Bullock 2018 Bloomington, Indiana
%
%% set up aparcAsegFile

if or(isstring(atlas),ischar(atlas))
    if ~exist(atlas,'file')
    %apaprently necessary for matlab?
    spaceChar={' '};
    %I dont know why this is necessary
    quoteString=strcat('mri_convert',spaceChar, fsDir,'/mri/',mgzfile ,spaceChar, niiPath);
    quoteString=quoteString{:};
    [status result] = system(quoteString, '-echo');
    if status~=0
        warning('/n Error generating aseg nifti file.  There may be a problem finding the file. Output: %s ',result)
        
    end
end
    atlas=niftiRead(atlas);
else
    %do nothing
end

%% run the merge roi function
[mergedROI] =bsc_ROIFromAtlasNums(atlas,ROInums, smoothKernel);

%% save it and ouput the nii and name

[ni, roiName]=dtiRoiNiftiFromMat (mergedROI,atlas,mergedROI.name,1);
end