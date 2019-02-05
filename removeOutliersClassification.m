function classification= removeOutliersClassification(classification,wbFG, centroidSD, lengthSD,maxIter,selectPrune)
%
% FiberIndexes= removeOutliersClassification(classification,wbFG, centroidSD, lengthSD, selectPrune)
%
% Performs outlier removal given input parameters. Defauts to 4 and 4 if
% nothing put in.
%
% INPUTS: 
%
% -classification: Either the path to structure or the structure itself.
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
% -wbFG: Either the path to the fe or WBFG, or the
%  structure itself
%
% -centroid SD: paramater for the sd deviance cut off from distance from
%  centroid of fibers.  Fibers outside of this distance will be pruned.
%
% -length SD: paramater for the sd deviance cut off for fiber length
%  fibers exceeding the threshold will be pruned.
%
% -maxIter: the maximum number of iterations to run.
%
% -selectPrune: a vector of indexes into the classification.names structure
%  indicating which tracts the user would like to have pruned
%% preliminaries

if ischar(classification)
    load(classification);
end

% loads file if a string was passed 
if ischar(wbFG)
    wbFG = load(wbFG);
    %if it is a fe structure, get the wbFG out of it
    if isfield(wbFG, 'fe')
        wbFG = feGet(wbFG.fe, 'fibers acpc');
    end
else
    if isfield(wbFG, 'fg')
        wbFG = feGet(wbFG, 'fibers acpc');
    end
end

%if it is a fe structure, get the wbFG out of it
if isfield(wbFG, 'fe')
    wbFG = feGet(wbFG.fe, 'fibers acpc');
end

% if selectPrune is not defined, just assume the user wants to prune all
% the tracts.
if notDefined('selectPrune')
selectPrune=1:length(classification.names);
end

if notDefined('centroidSD')
    centroidSD=4;
end

if notDefined('lengthSD')
    lengthSD=4;
end

%% begin outlier removal

%count number of classified tracts
classCount=length(find(classification.index>0));

%create blankFG
tractFG=wbFG;
tractFG.fibers=[];

pruneTotal=0;
for itracts=selectPrune
    tractFG.name=classification.names{itracts};
    indexes=find(classification.index==itracts);
    tractFG.fibers=wbFG.fibers(indexes);
    [~, keep]=mbaComputeFibersOutliers(tractFG,centroidSD,lengthSD,maxIter);
    
    %count number to be pruned
    trackPrune=length(indexes(~keep));   
    
    if trackPrune~=0
    fprintf('\n %i of %i streamlines (%.2f percent) pruned for the %s',trackPrune, length(keep), (trackPrune/length(keep))*100 ,tractFG.name)
    
    %find the indexes that were removed (i.e., indexes(~keep)) and set them
    %to zero, as though they had not been classified in the first place.
    classification.index(indexes(~keep))=0;
    
    %add to running total
    pruneTotal=(pruneTotal+trackPrune);
    else 
    fprintf('\n no streamlines removed for for the %s', tractFG.name)
    end
end

fprintf('\n %i total streamlines pruned (%4.2f percent of origional input)', pruneTotal, (pruneTotal/classCount)*100)

end
    
    
