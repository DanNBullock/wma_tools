function tractStruc = bsc_makeFGsFromClassification_v5(classification, wbfg)
%    tractStruc = bsc_makeFGsFromClassification(classification, wbfg)
%
%    Purpose:  This function creates a stucture containing all of the fiber
%    groups contained within a classification structure.  This facilitates
%    plotting or other visualization/analysis
%
%    INPUTS
%  -classification: Either the path to structure or the structure itself.
%   The strucure has a field "names" with (N) names of the tracts
%   classified while the field "indexes" has a j long vector (where  j =
%   the nubmer of streamlines in wbFG (i.e. length(wbFG.fibers)).  This j
%   long vector has a 0 for to indicate a streamline has gone unclassified,
%   or a number 1:N indicatate that the streamline has been classified as a
%   member of tract (N).
%
%   wbfg:  a structure containing the streamlines referenced in the
%   classification structure.  A whole brain fiber group, can be a tck
%   depending on your version of vistasoft.  Will load paths.
%
%
% OUTPUTS tractStruc:  structure wherein each row corresponds to a diffrent
% fg. Previous version added a layer of complexity to this, such there were
% fg and name fields.  This was redundant because there was a .name field
% within .fg.  As such, this is now just a cell structure (wherein each
% cell contains a fg object), and thus there are no fields.
%
%
% (C) Daniel Bullock, 2017, Indiana University
% Refactored 4/20/20
%% begin code

% loads object if path passed
if ischar(wbfg)
    wbfg = fgRead(wbfg);
else
    %do nothing
end

% loads classification structure if a path is passed
if ischar(classification)
    load(classification);
else
    %do nothing
end

% just for asthetic purposes, selects colors and applies to them to the
% tract
classificationGrouping = wma_classificationStrucGrouping(classification);
colorMapping = hsv(max(classificationGrouping.index));
% end

if isempty(classification.names)
    error('classification.names is empty.. something went wrong?');
end

% for each name in the classification.names structure finds the
% corresponding streamlines associated with that classification and creates
% an fg containg all relevant streamlines

%loops
for itracts=1:length(classification.names)
    %use dti toolset to make new holder with specified name
    tractStruc{itracts} = dtiNewFiberGroup(classification.names{itracts});
    %Set color, a bit more complicated than you'd expect
    currentColor=colorMapping(unique(classificationGrouping.index(classification.index==itracts)),:);
    tractStruc{itracts}.colorRgb=currentColor;
    %preface the creation of tract with feedback
    fprintf('\n creating tract for %s with %i streamlines', classification.names{itracts},sum(classification.index==itracts) );
    %assign streams to new fg object
    tractStruc{itracts}.fibers=wbfg.fibers(classification.index==itracts);
    %impose orientation on it
    tractStruc{itracts}=bsc_orientFgRAS_LPI(tractStruc{itracts});
end
end
