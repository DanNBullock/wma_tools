function tractStruc = bsc_makeFGsFromClassification_v2(classification, wbFG,coordScheme)
%    tractStruc = bsc_makeFGsFromClassification(classification, wbFG)
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
%   wbFG:  a structure containing the streamlines referenced in the
%   classification structure.  Can be either a fe structure or a whole
%   brain fiber group.  Will load paths.
%
%  coordScheme: the coordinate scheme you would like the output fibers in,
%  default is acpc.  Other option is IMG
%
% OUTPUTS tractStruc:  Cell array wherein each cell corresponds to a
% diffrent fg. Previous version added a layer of complexity to this.
%
% (C) Daniel Bullock, 2017, Indiana University
%
%% begin code

%probably won't work in long run.  Also not relevant in most cases.
%Apprently the previous assumption was that you'd only be using img space
%if you were inputting an fe structure.  
if notDefined('coordScheme')
    coordScheme='acpc';
end

% loads requisite structures from input
[wbFG, fe] = bsc_LoadAndParseFiberStructure(wbFG);

% loads classification structure if a path is passed
if ischar(classification)
    load(classification);
end

% if an fe structure is detected, alters classificaiton index to only
% include positively weighted fibers
if ~isempty(fe)
    classification=wma_clearNonvalidClassifications(classification,fe);
end

%find left and right tracts and assign them the same color.
classificationGrouping = wma_classificationStrucGrouping(classification);
colorMapping = distinguishable_colors(length(classificationGrouping.names),'k');

% for each name in the classification.names structure finds the
% corresponding streamlines associated with that classification and creates
% an fg containg all relevant streamlines
if strcmpi(coordScheme,'acpc')
    for itracts=1:length(classification.names)
        if ~isempty(unique(classificationGrouping.index(find(classification.index==itracts))));
            colorIndex=unique(classificationGrouping.index(find(classification.index==itracts)));
        else
            warning('\n Color assignment for tract %s skipped', classification.index{itracts});
        end
        
        tractStruc{itracts} = dtiNewFiberGroup(classification.names{itracts});
        tractStruc{itracts}.colorRgb=colorMapping(colorIndex,:);
        tractStruc{itracts}.fibers=wbFG.fibers(classification.index==itracts);
    end
elseif strcmpi(coordScheme,'img')
    for itracts=1:length(classification.names)
        
        if ~isempty(unique(classificationGrouping.index(find(classification.index==itracts))));
            colorIndex=unique(classificationGrouping.index(find(classification.index==itracts)));
        else
            warning('\n Color assignment for tract %s skipped', classification.index{itracts});
        end
        
        tractStruc{itracts} = dtiNewFiberGroup(classification.names{itracts});
        tractStruc{itracts}.colorRgb=colorMapping(colorIndex,:);
        tractStruc(itracts).fibers=fe.fg.fibers(classification.index==itracts);
    end
else
    fprintf('coordScheme input not understood')
end
end
