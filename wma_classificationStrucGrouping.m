function classificationGrouping = wma_classificationStrucGrouping(classification)




tractNames=classification.names;

%ostensible methods for label left and right
removalVec={'left','left ','l_','l ','right','right ','r_','r '};

flagMatrix=zeros(length(removalVec),length(tractNames));

%presumes that side indexing comes at front of name
for iRemoveLabels=1:length(removalVec)
    tempNames=cellfun(@(x) x(1:length(removalVec{iRemoveLabels})),tractNames,'UniformOutput',false);
    flagMatrix(iRemoveLabels,:)=~cellfun(@isempty,(strfind(lower(tempNames),removalVec{iRemoveLabels})));
end

%modify the flag vec for 'left' and 'left ' overlap
flagMatrix(1,:)=xor(flagMatrix(1,:),flagMatrix(2,:));
flagMatrix(5,:)=xor(flagMatrix(5,:),flagMatrix(6,:));

leftFlag=sum(flagMatrix(1:length(removalVec)/2,:),1)>0;
rightFlag=sum(flagMatrix((length(removalVec)/2+1):end,:),1)>0;

truncatedName=cell(1,length(tractNames));

for iNames=1:length(tractNames)
    for iRemoveLabels=1:length(removalVec)
        if flagMatrix(iRemoveLabels,iNames)
            
            truncatedName{iNames}=tractNames{iNames}(length(removalVec{iRemoveLabels})+1:end);
        end
    end
    if isempty(truncatedName{iNames})
        truncatedName{iNames}=tractNames{iNames};
    end
end

uniqueTruncated=unique(truncatedName,'stable');

classificationGrouping.names=uniqueTruncated;
classificationGrouping.index=zeros(length(classification.index),1);

for iUniqueNames=1:length(uniqueTruncated)
    bothIndexes=find(strcmp(uniqueTruncated{iUniqueNames},truncatedName));
    if length(bothIndexes)==2
        classificationGrouping.index(or(classification.index==bothIndexes(1),classification.index==bothIndexes(2)))=iUniqueNames;
    elseif length(bothIndexes)==1
        classificationGrouping.index(classification.index==bothIndexes)=iUniqueNames;
    else
        error('\n % i tracts are labeled as either right or left %s', length(bothIndexes),uniqueTruncated{iUniqueNames} );
        
    end
end
end
