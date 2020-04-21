function [outAtlas] = fnDeislandLabels_v3(atlas,maxisleSize,replaceVal)
% [olab] = fnDeislandLabels_v3(labels,outfile,maxisleSize,replaceVal)
%
%  This function uses bwareaopen to identify volumetric islands in
%  parcellations/atlases.  Useful for fixing freesurfer outputs, which are
%  known to have islands (i.e. 1 voxel, non-connected label instances).
%  NOTE: THE ISLAND VOXELS ARE REPALCED WITH THE replaceVal INPUT VARIABLE
%  Also note:  This function can be used to identify and replace problem
%  voxels.  After reassigning these voxels to replaceVal, you simply
%  selectively inflate into voxels identified as replaceVal.  See
%  bsc_inflateRelabelIslands for implementation.
%
%  INPUTS
%  atlas: either a path to a nfiti or the nifti object itself (a
%  parcellation or atlas, i.e. integer filled)
%
%  maxisleSize: the threshold size for islands in your atlas.  All islands
%  with fewer than this many voxels will be reassigned to replaceVal.
%
%  replaceVal:  islands will be reassigned this value.  Typically, input an
%  integer not already used in your atlas here
%
%  OUTPUTS
%  outAtlas: an modified version of the input atlas with the islands
%  removed
%
% voxles from outside of the largest continuous label (bad freesurfer labels)
%
%   Detailed explanation goes here
%
% Brent McPherson, (c) 2019, Indiana University
% Dan Bullock, slight modifications.
% Revised again 4/20/20

fprintf('fnDeislandLabels_v2 maxisleSize:%i replaceVal:%i\n', maxisleSize, replaceVal)

% read in aligned aparc+aseg .nii (read labels)
if or(ischar(atlas),isstr(atlas))
    atlas = niftiRead(atlas);
else
    %do nothing
end

% create a copy in output
outAtlas = atlas;

% get list of unique labels (0 is background)
uniqueLabels = unique(atlas.data(:));
uniqueLabels = uniqueLabels(uniqueLabels > 0);

% extract the data
data = atlas.data;
imgdim = size(atlas.data);

% preallocate output data labels
outData = zeros(imgdim);

%initialize number lost
nlost=0;
for iLabels = 1:size(uniqueLabels, 1)
    
    % create 3d binarized image of the label
    currentMask = int32(data == uniqueLabels(iLabels));
    
    % find islands of label THAT ARE LESS THAN THE MAX ISLAND SIZE
    % outputs a new binarized mask
    currentConnected = bwareaopen(currentMask,maxisleSize);
   
    % write the output labels from the largest component into the out data
    % object
    outData(currentConnected)=uniqueLabels(iLabels);
    
    %set the island voxels to have the replacement value
    outData(currentMask&~currentConnected)=replaceVal;
    nlost=nlost+length(find(currentMask&~currentConnected));
end
fprintf('%i total island voxels lost',nlost)
% assign output data
outAtlas.data = outData;

end
