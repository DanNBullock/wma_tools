function [indexBool] = bsc_extractStreamIndByName(classification,tractName)
% [indexBool] = bsc_extractStreamIndByName(classification,tractName)
%
% super straightforward function to extract boolean vector indexing
% relevant streams.  Probably not actually necessary, but needed in order
% to make roust against failure to find errors.

%%
%find index of name
nameIndex=strcmp(classification.names,tractName);

%set indexBool to 0, no index can be -1, so this should be false
indexBool=classification.index==-1;

if ~isempty(nameIndex)
    relevantIndexes=classification.index==nameIndex;
    indexBool(relevantIndexes)=true;
else
    %do nothing, index bool should remain all false
    warning('\n no streamlines found for %s found in classification structure',tractName)
end
end
    