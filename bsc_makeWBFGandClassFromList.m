function [mergedFG,classification]=bsc_makeWBFGandClassFromList(tckList)

classification.names=[];
classification.index=[];

% Amend name of tract in classification structure
for ii = 1:length(tckList)
    
    [~,fgName,~] = fileparts(tckList{ii}); 
    currFG=wma_loadTck(tckList{ii});
    if ii==1
        [mergedFG, classification]=bsc_mergeFGandClass(currFG,classification);
    else
        [mergedFG, classification]=bsc_mergeFGandClass(mergedFG,classification);
    end
    classification.names{ii} = fgName;
end