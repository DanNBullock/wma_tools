function [inflatedAtlas] =bsc_inflateLabels_v3(atlas,inflateItr)
% This function inflates the DK2009 atlas from freesurfer such that
% cortical labels are extended into the white matter and in to unknown
% labels
%
% Inputs:
% -atlas: path to THIS SUBJECT'S DK2009 Atlas in
% nii.gz format, or the object itself
% -inflateItr:  iterations of inflation.  Essentially akin to a radial
% smoothing kernel.

%
% OUTPUTS:
%  inflatedAtlas: the inflated version of the atlas.
% (C) Daniel Bullock, 2019, Indiana University
%
%% Begin code
disp('bsc_inflateLabels');

% loads object if path passed
if ischar(atlas)
atlas = niftiRead(atlas);
else
    %do nothing
end

%set relevant ROI indicies
greyMatterROIS=[[101:1:175]+12000 [101:1:175]+11000];
subcorticalROIS=[10:13 17:20 26 58 27 49:56 59 ];
%spineROIS=[16 28 60];
%cerebellumROIS=[8 47 7 46 ];
%ventricleROIS=[31 63 11 50 4 43 77 14 24 15 44 5 62 30 80];
wmROIS=[41 2];
%ccROIS=[251:255];
%unknownROIS=[0 2000 1000];
%OpticCROI=[85];


%set inflation targets
inflateTargROIs=[wmROIS 2000 1000 999];

%begin timer
%maybe a dupliate,but whatever
workingAtlas=atlas;
for iIter=1:inflateItr
    
    %inflate the grey and subcortical rois
    inflateROis=bsc_roiFromAtlasNums(workingAtlas,[greyMatterROIS subcorticalROIS], 3);
    noInflateROIS=bsc_roiFromAtlasNums(workingAtlas,[greyMatterROIS subcorticalROIS], 1);
    
    %create nifti of the inflated rois
    [inflateInfo, ~] = dtiRoiNiftiFromMat(inflateROis,atlas,'roisInflate',false);
    [noInflateInfo, ~] = dtiRoiNiftiFromMat(noInflateROIS,atlas,'roisNoInflate',false);
    
    %find the difference, this is the inflation margin
    differenceData=inflateInfo.data~=noInflateInfo.data;
    
    %go ahead and pad the data now
    paddedData=padarray(workingAtlas.data,[1 1 1]);
    %find the indexes that are both inflated into and targets.
    targetVoxelMask=differenceData & ismember(workingAtlas.data,inflateTargROIs);
    
    %for testing if you want, to see what is going to be replaced:
    %[counts, labels]=groupcounts(atlasIterCur.data(targetVoxelMask))
    
    %pad the mask you get from that
    paddedTargetVoxelMask=padarray(targetVoxelMask,[1 1 1]);
    %find the padded indexes
    targetVoxelPadedInd=find(paddedTargetVoxelMask);
    %get dimensions for the current image
    imgdim = size(paddedData);
    %taken from bsc_inflateRelabelIslands
    for iVoxels=1:length(targetVoxelPadedInd)
        currentVoxel=targetVoxelPadedInd(iVoxels);
        %convert to a coordinate
        [curX, curY, curZ]=ind2sub(imgdim,currentVoxel);
        %find the window of label values
        voxelVote=paddedData([curX-1:curX+1],[curY-1:curY+1],[curZ-1:curZ+1]);
        %count what groups there are
        [counts, labels]=groupcounts(voxelVote(~ismember(voxelVote,inflateTargROIs)));
        %determine vote winner
        voteWinner=labels(find(max(counts)));
        %change the output
        paddedData(currentVoxel)= voteWinner;
    end
    fprintf('(bsc_inflateLabels_v3.m) %i voxels for iteration %i\n',length(targetVoxelPadedInd),iIter)
    workingAtlas.data=paddedData(2:end-1,2:end-1,2:end-1);
end
%set the output
inflatedAtlas=workingAtlas;    
end  


