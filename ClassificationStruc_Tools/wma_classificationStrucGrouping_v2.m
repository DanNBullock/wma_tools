function classificationGrouping = wma_classificationStrucGrouping_v2(classification)
% classificationGrouping = wma_classificationStrucGrouping(classification)
%
% This function finds which pairs of names (and thus index labels)
% correspond to left/right variants of the same tract.  It presumes a
% finite number of possible lables that can indicate this status (as
% established in the removalVec variable.  It the creates another
% classificaiton structure, essentially unifying the tract into one, which
% has no left right designation.
%
%  Inputs
%  classification:  standardly organized classification
%
%  Outputs:
%  classificationGrouping:  a classification structure with the left and
%  right classifications merged
%
%  Dan Bullock 2019
%
%  Refactored in 4/2/2021 to incorporate regex.  Now gets around issues
%  relating to short names.  Also takes care of suffixes too.
%% begin code

%extract names vector
tractNames=classification.names;

%ostensible methods for label left and right
prefixVec={'left ','left_','left','l_','l ','right ','right_','right','r_','r '};

%create a blank structure to mark those structures which are changed
prefixFlagMatrix=zeros(length(prefixVec),length(tractNames));

%do the fix for prefixes
for iRemoveLabels=1:length(prefixVec)
    %extract current label
    currLabel=prefixVec{iRemoveLabels};
    
    %remove the prefix from each name, if it occurs at the front
    alteredNames=cellfun(@(x) regexprep(x,strcat('^',currLabel,'(?i)'),'','ignorecase'), tractNames,'UniformOutput', false);
    
    %check for differences
    newNames=setdiff(alteredNames,tractNames);
    
    %create blank vector to hold indexes of modified terms
    blankModVec=[];
    
    %loop across
    for iModifications =1:length(newNames)
        
        %find indexes of modified term
        modifiedIdexes=find(strcmp(newNames{iModifications},alteredNames));
        blankModVec=horzcat(blankModVec,modifiedIdexes);
    end
    
    %set matrix entries to true where relevant
    prefixFlagMatrix(iRemoveLabels,blankModVec)= true;
    
    %set names to new modified version
    tractNames=alteredNames;
end

%now do the same thing for suffixes
%ostensible methods for label left and right
suffixVec={' left','_left','left','_l',' l',' right','_right','right','_r',' r'};

%create a blank structure to mark those structures which are changed
suffixFlagMatrix=zeros(length(suffixVec),length(tractNames));

%do the fix for prefixes
for iRemoveLabels=1:length(suffixVec)
    %extract current label
    currLabel=suffixVec{iRemoveLabels};
    
    %remove the prefix from each name, if it occurs at the front
    alteredNames=cellfun(@(x) regexprep(x,strcat(currLabel,'$','(?i)'),'','ignorecase'), tractNames,'UniformOutput', false);
    
    %check for differences
    newNames=setdiff(alteredNames,tractNames);
    
    %create blank vector to hold indexes of modified terms
    blankModVec=[];
    
    %loop across
    for iModifications =1:length(newNames)
        
        %find indexes of modified term
        modifiedIdexes=find(strcmp(newNames{iModifications},alteredNames));
        blankModVec=horzcat(blankModVec,modifiedIdexes);
    end
    
    %set matrix entries to true where relevant
    suffixFlagMatrix(iRemoveLabels,blankModVec)= true;
    
    %set names to new modified version
    tractNames=alteredNames;
end

%modify the flag vec for 'left' and 'left ' overlap
%not sure what this is used for now
prefixFlagMatrix(1,:)=xor(prefixFlagMatrix(1,:),prefixFlagMatrix(2,:));
prefixFlagMatrix(5,:)=xor(prefixFlagMatrix(5,:),prefixFlagMatrix(6,:));

%get the unique names out
uniqueTruncated=unique(tractNames,'stable');

%create blank classification structure
classificationGrouping.names=uniqueTruncated;
classificationGrouping.index=zeros(length(classification.index),1);

%loop across the new unique names
for iUniqueNames=1:length(uniqueTruncated)
    %find the indexes 
    bothIndexes=find(strcmp(uniqueTruncated{iUniqueNames},tractNames));
    %if there are two, find the the .index entries that correspond to these
    %and merge them.
    if length(bothIndexes)==2
        classificationGrouping.index(or(classification.index==bothIndexes(1),classification.index==bothIndexes(2)))=iUniqueNames;
    %otherwise, not much else to do, just assign the indexes
    elseif length(bothIndexes)==1
        classificationGrouping.index(classification.index==bothIndexes)=iUniqueNames;
        %if there are more than 2, you have a problem
    else
        error('\n % i tracts are labeled as either right or left %s', length(bothIndexes),uniqueTruncated{iUniqueNames} ); 
    end
end

end