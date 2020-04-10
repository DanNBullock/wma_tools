function [classification] = bsc_spliceClassifications(baseClassification,classificationAdd)
% [classification] =bsc_reconcileClassifications(baseClassification,classificationAdd)
%
% DESCRIPTION: This function is for splicing two classification structures.  In
% splicing, the presumption is that these correspond to separate fg
% entities (i.e. collections of streamlines) but which may nonetheless
% correspond to the same anatomical and thus named structure
%
% INPUTS:
%
% baseClassification: a base classification structure
%
% classificationAdd: a classification structure that is to be added to the
%                    base classification structure.
%
% OUTPUTS:
%
% classification: the spliced / merged classificaiton 
%
%  (C) Daniel Bullock 2018 Bloomington
%% Begin code


baseNames=baseClassification.names;
baseNameNum=length(baseClassification.names);

addNames=classificationAdd.names;
addNameNum=length(classificationAdd.names);

presumeNameNum=baseNameNum+addNameNum;

uniqueNamesTotal=unique(horzcat(baseNames,addNames),'stable');
uniqueNamesLength=length(uniqueNamesTotal);

classification.names=[];
classification.index=zeros(length(baseClassification.index)+length(classificationAdd.index),1);

% if the name sets are exclusive, in that both classification structures
% feature entirely unique tracts
if uniqueNamesLength==presumeNameNum
    %smash the name vecs together
    classification.names=horzcat(baseNames,addNames);
    %adjust the indexing numerals of the classification index that is to be
    %added such that the indexes correspond to the new name vector.
    classificationAdd.index(classificationAdd.index>0)=classificationAdd.index(classificationAdd.index>0)+baseNameNum;
    
       
    %ok, but first actually put in the origional ones too
    classification.index(1:length(baseClassification.index))=baseClassification.index;
        
    %Basically, leave the dictionary values for the existing tracts (i.e. tracts
    %from baseClassification) alone, but change the dictionary values of the
    %streamline indexes corresponding to the classificationAdd input to
    %correspond to the newly modified values from classificationAdd.
    classification.index(find(classificationAdd.index>0)+length(baseClassification.index))=classificationAdd.index(classificationAdd.index>0);
else
    %in the event that there is some overap between the two classification
    %structures, such that some anatomical structures are found in both,
    %and some are found in only one, we now have to do things a bit
    %differently...
    
    classification.names=uniqueNamesTotal;
    
    for iNames=1:length(classification.names)
        baseNameIndex=find(strcmp(baseNames,uniqueNamesTotal{iNames}));
        addNameIndex=find(strcmp(addNames,uniqueNamesTotal{iNames}));
        
        if ~isempty(baseNameIndex)
            baseIndexes=find(baseClassification.index==baseNameIndex);
        else
            baseIndexes=[];
        end
        
        if ~isempty(addNameIndex)
            addIndexes=find(classificationAdd.index==addNameIndex);
        else
            addIndexes=[];
        end
        %for whatever reason, it seems like I understood the splicing issue
        %when I wrote this part, so no changes necessary here.
        classification.index(baseIndexes)=iNames;
        classification.index(addIndexes+length(baseClassification.index))=iNames;
        
    end
end
    
