function [booleanOut]=bsc_applyEndpointCriteria(fg, varargin)
% [booleanOut]=bsc_applyEndpointCriteria(fg, varargin)
% DESCRIPTION:
% This function evaluates a number of criteria relative to the endpoints of
% tracts and outputs a boolean vector indicating which endpoints meet all
% criteria.
%
% INPUTS:
% -fg: the fg the criteria will be applied to 
%
% -varargin: interprets sequential triplets of inputs such that
%    (1) a coordinate
%    (2) a relative criteria
%    (3) whether this criteria should apply to 'both', 'one', or 'neither' endpoint
%
%  ALTERNATIVE INPUT for #1
%
%  instead of a coordinate a plane can be input.  In such a case the plane
%  will be interpreted as the singular value of the plane (i.e. the "narrow"
%  coordinate).
%
% OUTPUTS:
% -booleanOut: boolean vector indicating which streamlines meet all criteria
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
        refCoord=varargin{iInputs*3-2};
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
        case 'either'
            currBool=or(sum(boolsOut,2)==1,sum(boolsOut,2)==2);
    end
    booleanOut=currBool&booleanOut;
    
end