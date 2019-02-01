function [booleanOut]=bsc_applyEndpointCriteria(fg, varargin)
% [booleanOut]=bsc_applyEndpointCriteria(fg, varargin)
% DESCRIPTION:
% This function evaluates a number of criteria relative to the endpoints of
% tracts and outputs a boolean vector indicating which endpoints meet all
% criteria.
%
% INPUTS:
% -fg:fg s Will extract
% midpoints from fg structure if that is put in instead.  It is recommended
% that midpoints be put in if this function is called multiple times within
% a script in order to increase speed.
%
% -varargin: interprets sequential pairings of inputs as (1) a coordinate
% and (2) a relative criteria.  For example
% bsc_applyEndpointCriteria(midpointsIn, [0 0 0], 'inferior') would output
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
allStreams=fg.fibers;
endpoints1=zeros(length(allStreams),3);
endpoints2=zeros(length(allStreams),3);
for iFibers=1:length(allStreams)

    %curStreamline=fg.fibers{iFibers};

    endpoints1(iFibers,:)=allStreams{iFibers}(:,1);
    endpoints2(iFibers,:)=allStreams{iFibers}(:,end);
end

booleanOut=true(1,length(fg.fibers))';

% redefine refCoord if plane
for iInputs=1:length(varargin)/3
    
    %if a plane is input, extract a single coordinate from it.
    if isstruct(varargin{iInputs*3-2})
        refCoord=[mean(varargin{iInputs*3-2}.coords(:,1)),mean(varargin{iInputs*3-2}.coords(:,2)),mean(varargin{iInputs*3-2}.coords(:,3))];
    else
        %no need to redfine refCoord
    end
    
    %determine if refcoord is on the left or the right.  Will cause problems
    %if it is right on the 0 line.
    LRflag=refCoord(:,1)<0;
    
    
    if any([strcmp(lower(varargin{iInputs*3-1}),'top'),...
            strcmp(lower(varargin{iInputs*3-1}),'superior')])
        dimension=3;
        valence=1;
    elseif any([strcmp(lower(varargin{iInputs*3-1}),'bottom'),...
            strcmp(lower(varargin{iInputs*3-1}),'inferior')])
        dimension=3;
        valence=-1;
    elseif strcmp(lower(varargin{iInputs*3-1}),'anterior')
        dimension=2;
        valence=1;
    elseif strcmp(lower(varargin{iInputs*3-1}),'posterior')
        dimension=2;
        valence=-1;
    elseif strcmp(lower(varargin{iInputs*3-1}),'medial')
        dimension=1;
        if LRflag
            valence=1;
        else
            valence=-1;
        end
    elseif strcmp(lower(varargin{iInputs*3-1}),'lateral')
        dimension=1;
        if LRflag
            valence=-1;
        else
            valence=1;
        end
    end
    
    criteriaCoords=[endpoints1(:,dimension) endpoints2(:,dimension)];
    
    if valence==1
        boolsOut=criteriaCoords>refCoord(dimension);
    else
        boolsOut=criteriaCoords<refCoord(dimension);
    end
    
    switch lower(varargin{iInputs*3})
        case 'one'
            currBool=sum(boolsOut,2)==1;
        case 'both'
            currBool=sum(boolsOut,2)==2;
        case 'neither'
            currBool=sum(boolsOut,2)==0;
    end
    booleanOut=currBool&booleanOut;
    
end