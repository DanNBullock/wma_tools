function [classificationOut] = bsc_splitTractAtPoint(fg, coordinate,location ,classification)
%[classificationOut] = bsc_splitTractAtPoint(fg, coordinate,location
%,classification)
%
% Description This function splits a tract at a given location based on the
% relative location of the closest node of a given streamline to a given
% plane.  As an illustration, this code is the generalized version of
% previous code which was used to separate the pArc and TPC from the same
% body at the bottom of the IPS.
%
% INPUTS
%
% fg:  either the fiber group corresponding to specific tract or (with the
% addition of a passed classification structure) a whole brain tractogram.
%
% coordinate: The XYZ coordinate of the point at which the split is to
% occur.
%
% location: specification of dimension division.  A string.  Either "x",
% "y", or "z"; or one of the anatomical relative position terms (see switch
% case - specificaiton of a single option, i.e. superior, necessarily
% implies a superior-inferior split) ;
%
%  ALTERNATE INPUTS
%
%  If a plane is put in for both coordinate and location then the following
%  will occur:
%
%  The plane input in coordinate will be used to extract the nearest node
%  from the relevant streamlines.  The streamlines will then be split
%  according to the location of those nodes relative to the plane input
%  into the location variable.  For example, a horizontal plane at the
%  bottom of the IPS would cause all streamlines to be judged at the node
%  that is closest to this height while a plane vertically and anteriorly
%  (and inferiorly and posteriorly) would cause the streamlines to be
%  devided in accordance with which sides of the plane the relevant nodes
%  were on.
%
%
% classification:  A single tract classification structure (i.e. from
% bsc_extractClassificatio).  If you pass a whole brain tractogram this
% will select
%
% OUTPUTS:
%
% classificationOut: a 2 tract classification structure.  Names inferred
% from relevant input.
%
% (C) Daniel Bullock 2018 Bloomington
%
%% Begin Code

if isstruct(coordinate) & isstruct(location)
    %if planes are passed it sets the variables accordingly. detect the
    %orientation of each
    uniqueCoordCountA=[length(unique(coordinate.coords(:,1))),length(unique(coordinate.coords(:,2))), length(unique(coordinate.coords(:,3)))];
    uniqueCoordCountB=[length(unique(location.coords(:,1))),length(unique(location.coords(:,2))), length(unique(location.coords(:,3)))];
    %with those counts determine the actual orientation of the plane
    narrowDimA=find(uniqueCoordCountA==min(uniqueCoordCountA));
    narrowDimB=find(uniqueCoordCountB==min(uniqueCoordCountB));
    %if planes are from same dimension throw error
    if narrowDimA==narrowDimB
        error('\n input planes are of same orientation')
    end
    %set planar coordinate for later node-find operation
    switch narrowDimA
        case 1
            refCoord=[mean(coordinate.coords(:,1)),0,0];
            planeDim=1;
        case 2
            refCoord=[0,mean(coordinate.coords(:,2)),0];
            planeDim=2;
        case 3
            refCoord=[0,0,mean(coordinate.coords(:,3))];
            planeDim=3;
    end
    
    %set planar coordinate for separation AND determine planar
    %separation
    switch narrowDimB
        case 1
            refCoord(1)=mean(location.coords(:,1));
            location='x';
            splitDim=1;
        case 2
            refCoord(2)=mean(location.coords(:,2));
            location='y';
            splitDim=2;
        case 3
            refCoord(3)=mean(location.coords(:,3));
            location='z';
            splitDim=3;
    end
    
elseif xor(isstruct(coordinate), isstruct(location))
    error('\n unrecognized input; mismatch between input types');
    
else
    refCoord=coordinate;
end

%if a classification structure is passed in
if exist('classification')
    fgExtract=bsc_makeFGsFromClassification_v2(classification, fg);
    fg=fgExtract{1};
else
    %do nothing
end
%determine if left or right point.  If tract is on left, x coords are less
%than zero, LR flag =1.
if refCoord(1)<0
    LRflag=1;
elseif refCoord(1)>0
    LRflag=0;
elseif refCoord(1)==0
    LRflag=NaN;
end

%holdover from older code, may be relevant later
switch lower(location)
    case 'top'
        splitDim=3;
    case 'superior'
        splitDim=3;
    case 'bottom'
        splitDim=3;
    case 'inferior'
        splitDim=3;
    case 'anterior'
        splitDim=2;
    case 'posterior'
        splitDim=2;
    case 'x'
        splitDim=1;
    case 'y'
        splitDim=2;
    case 'z'
        splitDim=3;
        %possible to deal with medial lateral issues here.
    case 'medial'
        if LRflag==1
            splitDim=1;
        elseif LRflag==0
            splitDim=1;
        elseif isnan(LRflag)
            splitDim=1;
        end
    case 'lateral'
        if LRflag==1
            splitDim=1;
        elseif LRflag==0
            splitDim=1;
        elseif isnan(LRflag)
            splitDim=1;
        end
end

%find node of interest
for iFibers = 1:length(fg.fibers)
    planeDist=abs(refCoord(planeDim) -fg.fibers {iFibers}(planeDim,:));
    nodeOfInterest=find(min(planeDist)==planeDist);
    %NOTE: the indexing (1) on nodeOfInterest picks the first node in the
    %case that two nodes are the same min distance.
    TopPoints(:,iFibers)=fg.fibers {iFibers}(:,nodeOfInterest(1));
end

%find fg indexes of relevant tracts
streamsGreater=find(TopPoints(splitDim,:)>refCoord(splitDim));
streamsLess=find(TopPoints(splitDim,:)<refCoord(splitDim));

%fill in new classificaiton structure
if exist('classification')
    classificationOut=classification;
else
    classificationOut.index=zeros(length(fg.fibers),1);
end

%assign the streams with the appropriate name given the dimension
switch splitDim
    case 1
        if exist('classification')
            wholeFGIndexes=find(classification.index);
            if LRflag==1
                %medial
                classificationOut.names{1}=strcat('medial_',fg.name);
                classificationOut.index(wholeFGIndexes(streamsGreater))=1;
                %lateral
                classificationOut.names{2}=strcat('lateral_',fg.name);
                classificationOut.index(wholeFGIndexes(streamsLess))=2;
            elseif LRflag==0
                %lateral
                classificationOut.names{1}=strcat('lateral_',fg.name);
                classificationOut.index(wholeFGIndexes(streamsGreater))=1;
                %medial
                classificationOut.names{2}=strcat('medial_',fg.name);
                classificationOut.index(wholeFGIndexes(streamsLess))=2;
            elseif isnan(LRflag)
                %Right
                classificationOut.names{1}=strcat('right_',fg.name);
                classificationOut.index(wholeFGIndexes(streamsGreater))=1;
                %Left
                classificationOut.names{2}=strcat('left_',fg.name);
                classificationOut.index(wholeFGIndexes(streamsLess))=2;
            end
        else
            if LRflag==1
                %medial
                classificationOut.names{1}=strcat('medial_',fg.name);
                classificationOut.index(streamsGreater)=1;
                %lateral
                classificationOut.names{2}=strcat('lateral_',fg.name);
                classificationOut.index(streamsLess)=2;
            elseif LRflag==0
                %lateral
                classificationOut.names{1}=strcat('lateral_',fg.name);
                classificationOut.index(streamsGreater)=1;
                %medial
                classificationOut.names{2}=strcat('medial_',fg.name);
                classificationOut.index(streamsLess)=2;
            elseif isnan(LRflag)
                %Right
                classificationOut.names{1}=strcat('right_',fg.name);
                classificationOut.index(streamsGreater)=1;
                %Left
                classificationOut.names{2}=strcat('left_',fg.name);
                classificationOut.index(streamsLess)=2;
            end
        end
    case 2
        if exist('classification')
            wholeFGIndexes=find(classification.index);
            %anterior
            classificationOut.names{1}=strcat('anterior_',fg.name);
            classificationOut.index(wholeFGIndexes(streamsGreater))=1;
            %posterior
            classificationOut.names{2}=strcat('posterior_',fg.name);
            classificationOut.index(wholeFGIndexes(streamsLess))=2;
        else
            %anterior
            classificationOut.names{1}=strcat('anterior_',fg.name);
            classificationOut.index(streamsGreater)=1;
            %posterior
            classificationOut.names{2}=strcat('posterior_',fg.name);
            classificationOut.index(streamsLess)=2;
        end
    case 3
        if exist('classification')
            wholeFGIndexes=find(classification.index);
            %superior
            classificationOut.names{1}=strcat('superior_',fg.name);
            classificationOut.index(wholeFGIndexes(streamsGreater))=1;
            %inferior
            classificationOut.names{2}=strcat('inferior_',fg.name);
            classificationOut.index(wholeFGIndexes(streamsLess))=2;
        else
            %anterior
            classificationOut.names{1}=strcat('superior_',fg.name);
            classificationOut.index(streamsGreater)=1;
            %posterior
            classificationOut.names{2}=strcat('inferior_',fg.name);
            classificationOut.index(streamsLess)=2;
        end
end
end
%end