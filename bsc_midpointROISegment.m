function  [midpointsBool] = bsc_midpointROISegment(wbfgORMidpoints,ROI)
%bsc_MidpointROISegment(wbfgORMidpoints,ROI)
%
% this function determines whether the midpoint of a set of streamlines
% intersects with an ROI
%
% Inputs:
% wbfgORMidpoints: either the midpoints of a fiber group or the fibergroup
% itself
% ROI: the roi to use
%
% Outputs:
% -midpointsBool:  a boolean vector indicating which streamlines meet the
% specified criteria
%
% (C) Daniel Bullock, 2019, Indiana University

distanceCriteria=.5;

if isstruct(wbfgORMidpoints)
    allStreams=wbfgORMidpoints.fibers;
    midpoints=zeros(length(allStreams),3);
    for iFibers=1:length(allStreams)
        fiberNodeNum=round(length(allStreams{iFibers})/2);
        curStreamline=allStreams{iFibers};
        midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    end
    % if a coordinate vector is input, make sure it is rotated correctly.
    % I.e. vertical column.
else
    midpointsInDim=size(wbfgORMidpoints);
    if ~midpointsInDim(2)==3
        midpoints=rot90(wbfgORMidpoints);
    else
       midpoints= wbfgORMidpoints;
    end
end

ROI=bsc_loadAndParseROI(ROI);

%TODO - handle a case when ROI.coords is empty (nearpoints will crash)
[~, distr1e1]=nearpoints(midpoints',ROI.coords');

midpointsBool=distr1e1<=distanceCriteria;
midpointsBool=midpointsBool';
end
