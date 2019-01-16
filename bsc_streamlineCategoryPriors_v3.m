function [classificationOut] =bsc_streamlineCategoryPriors_v3(wbfg, fsDir,inflateITer)
%
% [classificationOut] =bsc_segmentCorpusCallosum(wbfg, fsDir)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University

%% parameter note & initialization

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

classificationMid=classificationOut;

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

greyMatterROIS=[[101:1:175]+12000 [101:1:175]+11000];
subcorticalROIS=[10:13 17:20 26 58 27 49:56 59 ];
spineROIS=[16 28 60];
cerebellumROIS=[8 47 7 46 ];
ventricleROIS=[31 63 11 50 4 43 77 14 24 15 44 5 62 30 80];
wmROIS=[41 2];
ccROIS=[251:255];
unknownROIS=[0 2000 1000];
OpticCROI=[85];

interHemisphere=bsc_makePlanarROI(atlasPath,0, 'x');

for icategories=1:length(wbfg.fibers)
    curStream=wbfg.fibers{icategories};
endpoints1(:,icategories)=curStream(:,1);
endpoints2(:,icategories)=curStream(:,end);
end

[~, interHemiBool] = wma_SegmentFascicleFromConnectome(wbfg, [{interHemisphere}], {'and' }, 'dud');

if inflateITer>0
[inflatedAtlas] =bsc_inflateLabels(fsDir,inflateITer);
else
    inflatedAtlas=niftiRead(atlasPath);
end

[endpoints1Identity] =bsc_atlasROINumsFromCoords(inflatedAtlas,endpoints1,'acpc');
[endpoints2Identity] =bsc_atlasROINumsFromCoords(inflatedAtlas,endpoints2,'acpc');

termination1=[];
termination2=[];
%streamName{iStreams}=[];


for iStreams=1:length(wbfg.fibers)
    if     ~isempty(find(endpoints1Identity(iStreams)==greyMatterROIS))
        termination1{iStreams}='cortex';
    elseif ~isempty(find(endpoints1Identity(iStreams)==subcorticalROIS))
        termination1{iStreams}='subcortical';
    elseif ~isempty(find(endpoints1Identity(iStreams)==spineROIS))
        termination1{iStreams}='spinal';
    elseif ~isempty(find(endpoints1Identity(iStreams)==cerebellumROIS))
        termination1{iStreams}='cerebellum';
    elseif ~isempty(find(endpoints1Identity(iStreams)==ventricleROIS))
        termination1{iStreams}='ventricle';
    elseif ~isempty(find(endpoints1Identity(iStreams)==unknownROIS))
        termination1{iStreams}='unlabeled';
    elseif ~isempty(find(endpoints1Identity(iStreams)==wmROIS))
        termination1{iStreams}='whiteMatter';
    elseif ~isempty(find(endpoints1Identity(iStreams)==ccROIS))
        termination1{iStreams}='CorpusCallosum';
    elseif ~isempty(find(endpoints1Identity(iStreams)==OpticCROI))
        termination1{iStreams}='OpticChi';
    end
    
    if     ~isempty(find(endpoints2Identity(iStreams)==greyMatterROIS))
        termination2{iStreams}='cortex';
    elseif ~isempty(find(endpoints2Identity(iStreams)==subcorticalROIS))
        termination2{iStreams}='subcortical';
    elseif ~isempty(find(endpoints2Identity(iStreams)==spineROIS))
        termination2{iStreams}='spinal';
    elseif ~isempty(find(endpoints2Identity(iStreams)==cerebellumROIS))
        termination2{iStreams}='cerebellum';
    elseif ~isempty(find(endpoints2Identity(iStreams)==ventricleROIS))
        termination2{iStreams}='ventricle';
    elseif ~isempty(find(endpoints2Identity(iStreams)==unknownROIS))
        termination2{iStreams}='unlabeled';
    elseif ~isempty(find(endpoints2Identity(iStreams)==wmROIS))
        termination2{iStreams}='whiteMatter';
    elseif ~isempty(find(endpoints2Identity(iStreams)==ccROIS))
        termination2{iStreams}='CorpusCallosum';
    elseif ~isempty(find(endpoints2Identity(iStreams)==OpticCROI))
        termination2{iStreams}='OpticChi';
    end
    terminationNames=sort({termination1{iStreams} termination2{iStreams}});
    if interHemiBool(iStreams)
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2},'_interHemi');
    else
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2});
    end
end

uniqueNames=unique(streamName);

summarizeNames={'CorpusCallosum' 'unlabeled' 'OpticChi' 'ventricle' 'whiteMatter'};

% for icategories=1:length(uniqueNames)
% classificationMid=bsc_concatClassificationCriteria(classificationMid,uniqueNames{icategories},strcmp(streamName,uniqueNames{icategories}));
% end


for icategories=1:length(uniqueNames)
    if contains(uniqueNames{icategories},summarizeNames)
        for isummary=1:length(summarizeNames)
            summaryIndex(isummary)=contains(uniqueNames{icategories},summarizeNames{isummary});
        end
        summaryIndexSingle=find(summaryIndex);
        classificationOut=bsc_concatClassificationCriteria(classificationOut,summarizeNames{summaryIndexSingle(1)},contains(streamName,uniqueNames{icategories}));
        clear summaryIndex
    else
        classificationOut=bsc_concatClassificationCriteria(classificationOut,uniqueNames{icategories},contains(streamName,uniqueNames{icategories}));
    end
end

end