function [singleClassOUT]=bsc_extractClassification(classificationIN,nameOrIndex)
% [singleClassOUT]=bsc_extractClassification(classificationIN,nameOrIndex)
% DESCRIPTION:
% This function creates a new classification structure with the only tract
% being the specified tract or index.  The use case for this is to avoid
% having to pass around segmented tracts in order to modify them.  Instead
% you pass around the whole tractogram and the classification structure.
%
% INPUTS:
% -classificationIN: A classificaiton structure
%
% -nameOrIndex: either the name of the tract (as entered in the name field
% of the classification structure) or the index of that name (i.e. the
% index that corresponds to that name and the number that represets the
% tract in the .index field).
%
% OUTPUTS:
% -singleClassOUT: a single tract classification structure.  Names field
% will only contain the name corresponding to nameOrIndex input.  Index
% will only contain 1 (for identified streamlines) and 0 (for all other
% streamlines
%
%  (C) Daniel Bullock 2018 Bloomington
%% Begin Code
%
% switch case depending on if index or name is input
if isnumeric(nameOrIndex)
    %find name
    extractName=classificationIN.names{nameOrIndex};
    %find streamline indexes
    extractIndex=classificationIN.index==nameOrIndex;
elseif ischar(nameOrIndex)
    %get index from names field
    foundIndex=find(ismember(classificationIN.names,nameOrIndex));
    %do basically the same thing as above
    extractName=nameOrIndex;
    extractIndex=classificationIN.index==foundIndex;
end

%now make the new classification structure
singleClassOUT.names{1}=extractName;
singleClassOUT.index=extractIndex;

end


        