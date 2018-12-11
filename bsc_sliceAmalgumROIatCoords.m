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

if notDefined('space')
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

%check for monotonicity
monoTonicDim=or(sum(distances>0,2)==length(distances),sum(distances<0,2)==length(distances));
candidateDim=find(monoTonicDim);

if length(candidateDim>1)
    fprintf('\n More than one monotonic dimension found. Inferring segmentation dimension from maximal span covered.')
    DistancesSum=[sum(distances(1,:)),sum(distances(2,:)),sum(distances(3,:))]';
    segDim=find(DistancesSum==max(DistancesSum(candidateDim)));
elseif length(candidateDim==1)
    segDim=candidateDim;
elseif length(candidateDim==0)
    error('no clear monotonic direction in input coordinate scheme.  Please consider reordering input coordinates.')
end

%find the dimension along which the coordinates are arranged

% assumes they are sequential in order to segment

%get the ROIs that are associated with the coordinates


%check for unknown or white matter.  IT WOULD SEEM THAT IN THE CURRENT
%FORM, THE MAT APPLICATION DOES NOT PROVIDE LABELS FOR EMPTY SPACE OR FOR
%WHITE MATTER.  THIS LEADS TO ALL NON CORTICAL AREAS BEING LABELED 0.
%GIVEN THAT WE DO NOT WANT TO CREATE GIANT ROIS THAT COVER THE ENTIRE
%IMAGE, WE WILL NEED TO EXCLUDE THEM.  Conceivably one solution to this
%would be to determine closest label from point. Will consider this update
%later.

[ROInums] =bsc_atlasROINumsFromCoords(atlasNifti,coords,space);
if ~ isempty(find(ROInums==0))
    troubleCoords=[coords(:,find(ROInums==0))']
    warning('coordinate %i has returned an unlabeled index.  Excluding from amalgum creation.',find(ROInums==0))
    fprintf('\nRemaining ROIs will still be sliced with ALL coordinates (including trouble coordinate)')
    ROInums=ROInums(~ROInums==0);
end
    
%obtain merged ROI from the obtained numbers
[mergedROI]  =bsc_roiFromAtlasNums(atlasNifti,ROInums, []);

%set key anatomical relations for subsegmentation
if segDim==1
    if coords(segDim,1)>coords(segDim,end) & coords(segDim,1)>0
        firstWord='lateral';
        lastWord='medial';
    else
        firstWord='medial';
        lastWord='lateral';
    end
    if coords(segDim,1)>coords(segDim,end)& coords(segDim,1)<0
        firstWord='medial';
        lastWord='lateral';
    else
        firstWord='lateral';
        lastWord='medial';
    end
elseif segDim==2
    if coords(segDim,1)>coords(segDim,end)
        firstWord='anterior';
        lastWord='posterior';
    else
        firstWord='posterior';
        lastWord='anterior';
    end
elseif segDim==3
    if coords(segDim,1)>coords(segDim,end)
        firstWord='superior';
        lastWord='inferior';
    else
        firstWord='inferior';
        lastWord='superior';
    end
end

% set midpoint and total size
totalSlice=chkSize(2)+1;
HalfPoint=(totalSlice/2);

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
    if (ceil(HalfPoint)==icoords & ~(HalfPoint==icoords))
        mergedROI.name
        ROIs{icoords}.name=strcat(mergedROI.name,'_','middle');
    elseif icoords<=HalfPoint
        ROIs{icoords}.name=strcat(mergedROI.name,'_',firstWord,num2str(icoords));
    elseif icoords>HalfPoint
        ROIs{icoords}.name=strcat(mergedROI.name,'_',lastWord,num2str(icoords));
        
    end

    %amalgamate the remainder ROI
    [completedROI] = bsc_mergeROIs(ROIs{icoords}, completedROI);

end

[lastROI] = bsc_subtractROIs(completedROI, mergedROI);
ROIs{icoords+1}=lastROI;
ROIs{icoords+1}.name=strcat(mergedROI.name,'_',lastWord,num2str(icoords+1));

end
