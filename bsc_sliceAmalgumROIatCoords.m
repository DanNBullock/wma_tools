function [ROIs] =bsc_sliceAmalgumROIatCoords(atlasNifti,coords,space)
% [ROIs] =bsc_sliceAmalgumROIatCoords(atlasNifti,coords,space)
%
%  Purpose:  given a series of coordinates, find the atlas ROIs that
%  correspond to them, merge them in to one omnibus ROI, and then
%  sequentially slice it into regions defined by the coordinates.
%  Autodetects the orientation of the corrdinate sequence.
%
%  INPUTS:
%  -atlasNifti:  path to an atlas nifti.  An object works too.
%
%  -coords:  a 3 by N (where N is the number of coordinates) vector that
%  you would like to know the atlas ROI numbers for.
%
%  -space:  the space that the coordinates are in, either 'img' or 'acpc'.
%
%  OUTPUTS:
%
%  ROIs: an N +1 long cell structure with ROIs corresponding to 
%
% % (C) Daniel Bullock 2018 Bloomington, Indiana
%% Preliminaries 

if isempty(space)
    space='acpc'
end

%load atlas if necessary
if or(isstring(atlasNifti),ischar(atlasNifti))
    atlasNifti=niftiRead(atlasNifti);
else
    %do nothing
end

%reorient coords if necessary
chkSize=size(coords);

if ~chkSize(1)==3
    coords=rot90(coords);
else
    %do nothing
end


%% detect orientation of coordinates
chkSize=size(coords);
for icoords=1:chkSize(2)-1
    distances(1,icoords)=coords(1,icoords)-coords(1,icoords+1);
    distances(2,icoords)=coords(2,icoords)-coords(2,icoords+1);
    distances(3,icoords)=coords(3,icoords)-coords(3,icoords+1);
end

%find the dimension along which the coordinates are arranged
DistancesSum=[sum(distances(1,:)),sum(distances(2,:)),sum(distances(3,:))]';
maxDim=find(DistancesSum==max(DistancesSum));

% assumes they are sequential in order to segment

%get the ROIs that are associated with the coordinates
[ROInums] =bsc_atlasROINumsFromCoords(atlasNifti,coords,space);

%obtain merged ROI from the obtained numbers
[mergedROI]  =bsc_roiFromAtlasNums(atlasNifti,ROInums, []);

%set key anatomical relations for subsegmentation
if maxDim==1
    firstWord='lateral';
    lastWord='medial';
elseif maxDim==2
    firstWord='anterior';
    lastWord='posterior';
elseif maxDim==3
    firstWord='superior';
    lastWord='inferior';
end

% set midpoint and total size
totalSlice=chkSize(2)+1;
HalfPoint=(totalSlice/2)+.5;

%initiate roi to keep track of what has been completed (i.e. the 'other'
%boundary of the operation, other than the coordinate)
completedROI=mergedROI;
completedROI.coords=[];

%iterate slicing
for icoords=1:chkSize(2)
    %subtract what has been done from the amalgum ROI
    [remainderROI] = bsc_subtractROIs(completedROI, mergedROI);
    
    %Obtain the next slice
    ROIs{icoords}=bsc_modifyROI(atlasNifti,remainderROI, coords(:,icoords), firstWord);
    
    %name ROI appropriately
    if icoords<HalfPoint
        ROIs{icoords}.name=strcat(firstWord,num2str(icoords));
    elseif icoords>HalfPoint
        ROIs{icoords}.name=strcat(lastWord,num2str(icoords));
    elseif icoords>HalfPoint
        ROIs{icoords}.name='middle';
    end

    %amalgamate the remainder ROI
    [completedROI] = bsc_mergeROIs(ROIs{icoords}, completedROI);
    
end
