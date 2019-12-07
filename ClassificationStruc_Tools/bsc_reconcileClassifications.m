function [classification] = bsc_reconcileClassifications(baseClassification,classificationAdd)
%[classification] =bsc_reconcileClassifications(baseClassification,classificationAdd)
%
% DESCRIPTION: This function is for reconciling two classification
% structures.  In reconciling, the presumption is that these correspond to
% the *same* fg entities (i.e. collections of streamlines / source
% tractogram) but entities classified may be the same or different in the
% two classification structures.
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
% classification: the reconciled / merged classificaiton 
%
%  (C) Daniel Bullock 2018 Bloomington
%% Begin Code


if ~length(baseClassification.index)==length(classificationAdd.index)
    warning('\n Classification reconciliation assumption violated: input classifications of different lengths')
end
baseNames=baseClassification.names;
baseNameNum=length(baseClassification.names);

addNames=classificationAdd.names;
addNameNum=length(classificationAdd.names);

presumeNameNum=baseNameNum+addNameNum;

uniqueNamesTotal=unique(horzcat(baseNames,addNames),'stable');
uniqueNamesLength=length(uniqueNamesTotal);

classification.names=[];
%are there problems with setting the new clasification just equal to the
%base?
classification.index=baseClassification.index;


if uniqueNamesLength==presumeNameNum
    %hyper inelegant, guarenteed to cause problems.
    classification.names=horzcat(baseNames,addNames);
    classificationAdd.index(classificationAdd.index>0)=classificationAdd.index(classificationAdd.index>0)+baseNameNum;
    
    
    classification.index(classificationAdd.index>0)=classificationAdd.index(classificationAdd.index>0);
    
else
    
    classification.names=uniqueNamesTotal;
    for iNames=1:length(classification.names)
        baseNameIndex=find(strcmp(baseNames,uniqueNamesTotal{iNames}));
        addNameIndex=find(strcmp(addNames,uniqueNamesTotal{iNames}));
        
        if ~isempty(baseNameIndex)
            baseIndexes=find(baseClassification.index==baseNameIndex);
            unique(baseClassification.index);
        else
            baseIndexes=[];
        end
        
        if ~isempty(addNameIndex)
            addIndexes=find(classificationAdd.index==addNameIndex);
            unique(classificationAdd.index);
        else
            addIndexes=[];
        end
        
        %same fg structure assumed, so no need to do the splicing.
        classification.index(baseIndexes)=iNames;
        classification.index(addIndexes)=iNames;
        unique(classification.index);
    end
end
    