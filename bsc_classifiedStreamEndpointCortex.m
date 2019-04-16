function bsc_classifiedStreamEndpointCortex(wbFG, classification, fsDir, saveDir,subSelect, decayFunc, decayRadiusThresh)


%   bsc_plotClassifiedStreams(wbFG, classification, t1, view, saveDir,subSelect,colors)
%
%   PURPOSE: This function plots classified fibers using
%   mbaDisplayConnectome.  Will either plot all classified fiber tracts or
%   those that are subselected. Also, if you pass in an fe structure it
%   will only plot the validated fibers.  If you'd like to prune your
%   fibers as well, feel free to use removeOutliersClassification
%
%  -wbFG:  a structure containing the streamlines referenced in the
%    classification structure.  Can be either a fe structure or a whole
%    brain fiber group.  Will load paths.  SHOULD BE IN IMAGE SPACE.
%
%  -classification: Either the path to structure or the structure itself.
%   The strucure has a field "names" with (N) names of the tracts classified
%   while the field "indexes" has a j long vector (where  j = the nubmer of
%   streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%   a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%   indicatate that the streamline has been classified as a member of tract
%   (N).
%
%  -saveDir:  the directory you would like these figures saved to.  If not
%   defined, saves to current directory.
%
%  -subSelect: a vector corresponding to the indexes of the tracts (in the
%   classification.names structure) which you would like to plot.  If this
%   is not defined, then the function will plot all classified fiber
%   tracts.
%
%   -decayFunc: the decay function to use to calculate weight to throw in
%   voxels
%          ->uniform: no distance loss until threshold, then zero
%          ->linear: normalized linear cost to distance, 0 at dist > threshold
%          ->exponential: exponential cost to distance, 0 at dist > threshold
%          ->exact:  only counts exact endpoint voxel, uses floor()
%
%  -decayRadiusThresh: A distance threshold between tract endpoint and gray
%   matter voxels (default: 3 mm).
% 
%
% (C) Daniel Bullock, 2017, Indiana University


%% preliminaries
% loads requisite structures from input
[wbFG, fe] = bsc_LoadAndParseFiberStructure(wbFG);

%loads classificaiton file if a path is passed
if ischar(classification)
    load(classification);
end

fprintf('\n A total of %i streamlines are classified in this structure, representing %2.2f of the identified tracts ', length(find(classification.index>0)), length(find(classification.index>0))/length(find(fe.life.fit.weights))*100)

% if an fe structure is detected, alters classificaiton index to only
% include positively weighted fibers
if ~isempty(fe)
classification=wma_clearNonvalidClassifications(classification,fe);
end
fprintf('\n A total of %i streamlines are classified in this structure, representing %2.2f of the identified tracts ', length(find(classification.index>0)), length(find(classification.index>0))/length(find(fe.life.fit.weights))*100)

% if an fe structure is detected, alters classificaiton index to only
% include positively weighted fibers

fprintf('\n Of those, %i streamlines also have evidence, representing %2.2f of the streamlines with evidence.', length(find(classification.index>0)), length(find(classification.index>0))/length(find(fe.life.fit.weights>0))*100)
if ~isempty(fe)
classification=wma_clearNonvalidClassifications(classification,fe);
end

fprintf('\n Of those, %i streamlines also have evidence, representing %2.2f of the streamlines with evidence.', length(find(classification.index>0)), length(find(classification.index>0))/length(find(fe.life.fit.weights>0))*100)


% if user does not pass in a subselection
if notDefined('subSelect')
    subSelect=1:length(classification.names);
end

if notDefined('decayFunc'), decayFunc='box';end

if notDefined('decayRadiusThresh'), decayRadiusThresh='3';end

fprintf('\n using %s smoothing kernel with radius %i', decayFunc,decayRadiusThresh);


%% endpointMapping preliminaries


% Default is to compute sum of all endpoints within threshold distance


for iFGs = 1:length(subSelect)
    if ~isempty(find(classification.index==subSelect(iFGs)))
    %not really catching any of this, as it is being saved down
    [nii, nii2]=wma_endpointMapsDecay_v5(subSelect(iFGs), classification,fe,wbFG,fsDir, decayRadiusThresh, decayFunc);
    
    niftiWrite(nii,fullfile(saveDir,nii.fname));
    niftiWrite(nii2,fullfile(saveDir,nii2.fname));
    end
end

%deletes pool and temporary directory if it exists.



end