function [booleanOut]=bsc_applyMidpointCriteria(midpointsIn, varargin)
% [booleanOut]=bsc_applyMidpointCriteria(midpointsIn, varargin)
% DESCRIPTION:
% This function evaluates a number of criteria relative to the midpoints of
% tracts and outputs a boolean vector indicating which midpoints meet all
% criteria.
%
% INPUTS:
% -midpointsIn: path to THIS SUBJECT'S freesurfer directory.  Will extract
% midpoints from fg structure if that is put in instead.  It is recommended
% that midpoints be put in if this function is called multiple times within
% a script in order to increase speed.
%
% -varargin: interprets sequential pairings of inputs as (1) a coordinate
% and (2) a relative criteria.  For example
% bsc_applyMidpointCriteria(midpointsIn, [0 0 0], 'inferior') would output
% true values for all midpoints below y=0 (in acpc space).
%
%  ALTERNATIVE INPUT
%
%  instead of a coordinate a plane can be input.  In such a case the plane
%  will be interpreted as the singular value of he plane (i.e. the "narrow"
%  coordinate).
%
% OUTPUTS:
% -booleanOut: boolean vector indicating which midpoints meet all criteria
%
%  (C) Daniel Bullock 2018 Bloomington
%
%%  Begin Code
% if an fg structure is input, extract the midpoints
if isstruct(midpointsIn)
    allStreams=midpointsIn.fibers;
    midpoints=zeros(length(allStreams),3);
    for iFibers=1:length(allStreams)
        fiberNodeNum=floor(length(allStreams{iFibers})/2);
        curStreamline=allStreams{iFibers};
        midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    end
    % if a coordinate vectore is input, make sure it is rotated correctly.
    % I.e. vertical column.
else
    midpointsInDim=size(midpointsIn);
    if ~midpointsInDim(2)==3;
        midpoints=rot90(midpointsIn);
    else
       midpoints= midpointsIn;
    end
end

% redefine refCoord if plane
for iInputs=1:length(varargin)/2

%if a plane is input, extract a single coordinate from it.
if isstruct(varargin{iInputs*2-1})
    refCoord=[mean(varargin{iInputs*2-1}.coords(:,1)),mean(varargin{iInputs*2-1}.coords(:,2)),mean(varargin{iInputs*2-1}.coords(:,3))];
else
    %no need to redfine refCoord
end

%determine if refcoord is on the left or the right.  Will cause problems
%if it is right on the 0 line.
LRflag=refCoord(:,1)<0;

%execute switch statement and find coords that meet criteria
switch lower(varargin{iInputs*2})
    case 'top'
        refCoordSingle=refCoord(3);
        boolVec=midpoints(:,3)>refCoordSingle;
    case 'superior'
        refCoordSingle=refCoord(3);
        boolVec=midpoints(:,3)>refCoordSingle;
    case 'bottom'
        refCoordSingle=refCoord(3);
        boolVec=midpoints(:,3)<refCoordSingle;
    case 'inferior'
        refCoordSingle=refCoord(3);
        boolVec=midpoints(:,3)<refCoordSingle;
    case 'anterior'
        refCoordSingle=refCoord(2);
        boolVec=midpoints(:,2)>refCoordSingle;
    case 'posterior'
        refCoordSingle=refCoord(2);
        boolVec=midpoints(:,2)<refCoordSingle;
    case 'medial'
        if LRflag
            refCoordSingle=refCoord(1);
            boolVec=midpoints(:,1)>refCoordSingle;
        else
            refCoordSingle=refCoord(1);
            boolVec=midpoints(:,1)<refCoordSingle;
        end
    case 'lateral'
        if LRflag
            refCoordSingle=refCoord(1);
            boolVec=midpoints(:,1)<refCoordSingle;
        else
            refCoordSingle=refCoord(1);
            boolVec=midpoints(:,1)>refCoordSingle;
        end
        
end

% concat boolean vectors 
if iInputs==1
    booleanOut=boolVec;
else
booleanOut=and(booleanOut,boolVec);
end

end
