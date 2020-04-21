function outAtlas=bsc_inflateRelabelIslands(atlas)
%  fixedAtlas=bsc_inflateRelabelIslands
%
%  This function uses Brent Mcpherson's fnDeislandLabels to replace label
%  islands in an atlas/parcellation.  Islands are small, disconnected sets
%  of voxels that are not contiguous with the  major body of a label.  This
%  can cause probems if you are attempting to determine the borders of an
%  roi in an atlas.
%
%  INPUTS
%  atlas: either a path to a nfiti or the nifti object itself (a
%  parcellation or atlas, i.e. integer filled)
%
%  OUTPUTS
%  outAtlas:  An atlas with (hopefully) all of its islands removed
%  
%  Some rights reserved, Dan Bullock 4/20/20
%% begin code

if isstr(atlas)
    atlas = niftiRead(atlas);
else
    %do nothing
end

%change the label of islands with less than 5 voxels to 999
[outAtlas] = fnDeislandLabels_v3(atlas,7,999);

%from bsc_inflateLabels, except we are only doing the problem labels, 999.
%Theoretically you could also add other voxel labels in here (i.e. 2000
%1000 for DK2009) if you wanted to include them in the replacement process
trashLabel=[999];

%dangerous while loop, continues so long as a 999 is detected next to a
%labeled voxel

%get a mask of the island category
islandMask=outAtlas.data==999;
%find out how many individual islands there are
islandList= bwconncomp(islandMask);

%get dimensions of image
imgdim = size(atlas.data);

%use voting method to fix islands? iterating over islands might save time.
%AFTER TESTING APPEARS TO BE LIGNTHING FAST.
for iIslands=1:islandList.NumObjects
    % create a blank object for this island
    currentVoxels=islandList.PixelIdxList{iIslands};
    %iterate over the voxels
    voteVec=[];
    for iVoxels=1:length(currentVoxels)
        %get the ind for the current voxel
       currentVoxel=currentVoxels(iVoxels);
       %convert to a coordinate
       [curX, curY, curZ]=ind2sub(imgdim,currentVoxel);
       %probably shouldn't need to worry about padding
       voxelVote=outAtlas.data([curX-1:curX+1],[curY-1:curY+1],[curZ-1:curZ+1]);
    %concatonate the voxel votes  
    voteVec=vertcat(voteVec,voxelVote(:));
    end
    %get the counts for each label
    [counts, labels]=groupcounts(voteVec(all(voteVec~=[trashLabel 0],2)));
    %the one iwth the most votes wins
    voteWinner=labels(find(max(counts)));
    %change the output
    outAtlas.data(currentVoxels)=voteWinner;
end
end