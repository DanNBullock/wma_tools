function [classificationOUT]=bsc_concatClassificationCriteria(classificationIN,name,varargin)
% [classification]=bsc_concatClassificationCriteria(classification,name,varargin)
% DESCRIPTION:
% This function applies a series of AND opperations to the boolean vectors
% entered in to varargin, and then adds a tract classification to the
% passed classification structure, using the input name is the relevant
% tract name.  The presumption/use case for this is that the input boolean
% vectors represent a number of criteria that have been applied to a source
% tractogram such that the TRUE values are associated with streamlines that
% meet the criteria and FALSE values do not.  In this way, the resultant
% tract classification added to the classification structure represents
% those streamlines which meet all criteria.
%
% INPUTS:
% -classificationIN: A classificaiton structure
%
% -name: the name of the tract that the multiple criteria correspond to
%
% -varargin= a series of boolean vectors corresponding to various criteria
% checks (i.e. output from wma_SegmentFascicleFromConnectome,
% bsc_tractByEndpointROIs, etc.
%
%
% OUTPUTS:
% -classificationOUT: the updated classification structure
%
%  (C) Daniel Bullock 2018 Bloomington
%% Begin Code

%set initial criteria
omnibusCriteria=varargin{1};

%pull orientatiion from first input to impose on subsequent
orientationSet=find(size(omnibusCriteria)==max(size(omnibusCriteria)));

%set it now so that it will comply with classification orientation
if ~orientationSet==1
    omnibusCriteria=rot90(omnibusCriteria);
    orientationSet=1;
end

if length(varargin)>1
    for iCriteria=2:length(varargin)
        %determine current criteria
        curCriteriaOrientation=find(size(varargin{iCriteria})==max(size(varargin{iCriteria})));
        %rotate if necessary
        if curCriteriaOrientation~=orientationSet
            varargin{iCriteria}=rot90(varargin{iCriteria});
        else
            %no action necessary if not
        end
        omnibusCriteria=and(omnibusCriteria,varargin{iCriteria});
    end
else
    %no conjunction necessary if only one criteria passed in
end

%determine streamline number of source tractogram
if ~isempty(classificationIN)
    sourceLength=length(classificationIN.index);
    
    %create new classification structure just for this tract
    newClassification=[];
    newClassification.names{1}=name;
    newClassification.index=zeros(sourceLength,1);
    newClassification.index(omnibusCriteria)=true;
    
    %run the reconciliation
    classificationOUT=bsc_reconcileClassifications(classificationIN,newClassification);
    
else
    %if there is no input classification structure, just infer its
    %properties from the input booleans.
    sourceLength=length(omnibusCriteria);
    classificationOUT=[];
    classificationOUT.names{1}=name;
    classificationOUT.index=zeros(sourceLength,1);
    classificationOUT.index(omnibusCriteria)=true;
end

end