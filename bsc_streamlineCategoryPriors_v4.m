function [classificationOut] =bsc_streamlineCategoryPriors_v4(wbfg, fsDir,inflateITer)
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
fprintf('\n creating categorical segmentation')
wbfg
allStreams=wbfg.fibers;
%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(allStreams),1);
classificationOut

classificationMid=classificationOut;

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz')

greyMatterROIS=[[101:1:175]+12000 [101:1:175]+11000];
subcorticalROIS=[10:13 17:20 26 58 27 49:56 59 ];
spineROIS=[16 28 60];
cerebellumROIS=[8 47 7 46 ];
ventricleROIS=[31 63 11 50 4 43 77 14 24 15 44 5 62 30 80];
wmROIS=[41 2];
ccROIS=[251:255];
unknownROIS=[0 2000 1000];
OpticCROI=[85];

FrontalROIs=[[124 148 118 165 101 154 105 115 154 155 115 170 129 146 153 ...
    164 106 116 108 131 171 112 150 104 169 114 113 116 107 163 139 132 140]+11000 [124 148 118 165 101 154 105 115 154 155 115 170 129 146 153 ...
    164 106 116 108 131 171 112 150 104 169 114 113 116 107 163 139 132 140]+12000] ;

TemporalROIs=[[144 134 138 137 173 174 135 175 121 151 123 162 133]+11000 [144 134 138 137 173 174 135 175 121 151 123 162 133]+12000];

OccipitalROI=[[120 119 111 158 166 143 145 159 152 122 162 161 121 160 102]+11000 [120 119 111 158 166 143 145 159 152 122 162 161 121 160 102]+12000];

ParietalROI=[[157 127 168 136 126 125 156 128 141 172 147 109 103 130 110]+11000 [157 127 168 136 126 125 156 128 141 172 147 109 103 130 110]+12000];

pericROI=[[167]+11000 [167]+12000];

insulaROI=[[117 149]+11000 [117 149]+12000];
    
    
interHemisphere=bsc_makePlanarROI(atlasPath,0, 'x');


for icategories=1:length(allStreams)
    curStream=allStreams{icategories};
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


for iStreams=1:length(allStreams)
    if     ~isempty(find(endpoints1Identity(iStreams)==FrontalROIs))
        termination1{iStreams}='frontal';
    elseif ~isempty(find(endpoints1Identity(iStreams)==TemporalROIs))
        termination1{iStreams}='temporal';
    elseif ~isempty(find(endpoints1Identity(iStreams)==OccipitalROI))
        termination1{iStreams}='occipital';
    elseif ~isempty(find(endpoints1Identity(iStreams)==ParietalROI))
        termination1{iStreams}='parietal';
    elseif ~isempty(find(endpoints1Identity(iStreams)==pericROI))
        termination1{iStreams}='pericallosal';
    elseif ~isempty(find(endpoints1Identity(iStreams)==subcorticalROIS))
        termination1{iStreams}='subcortical';
    elseif ~isempty(find(endpoints1Identity(iStreams)==spineROIS))
        termination1{iStreams}='spinal';
    elseif ~isempty(find(endpoints1Identity(iStreams)==insulaROI))
        termination1{iStreams}='insula';
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
    
    if     ~isempty(find(endpoints2Identity(iStreams)==FrontalROIs))
        termination2{iStreams}='frontal';
    elseif ~isempty(find(endpoints2Identity(iStreams)==TemporalROIs))
        termination2{iStreams}='temporal';
    elseif ~isempty(find(endpoints2Identity(iStreams)==OccipitalROI))
        termination2{iStreams}='occipital';
    elseif ~isempty(find(endpoints2Identity(iStreams)==ParietalROI))
        termination2{iStreams}='parietal';
    elseif ~isempty(find(endpoints2Identity(iStreams)==pericROI))
        termination2{iStreams}='pericallosal';
    elseif ~isempty(find(endpoints2Identity(iStreams)==subcorticalROIS))
        termination2{iStreams}='subcortical';
    elseif ~isempty(find(endpoints2Identity(iStreams)==spineROIS))
        termination2{iStreams}='spinal';
    elseif ~isempty(find(endpoints2Identity(iStreams)==insulaROI))
        termination2{iStreams}='insula';
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
    if length(termination1)==length(termination2)
        terminationNames=sort({termination1{iStreams} termination2{iStreams}});
    else
        endpoints1Identity(iStreams)
        endpoints2Identity(iStreams)
        error('streamline identity unaccounted for')        
    end
    if interHemiBool(iStreams)
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2},'_interHemi');
    else
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2});
    end
end

uniqueNames=unique(streamName);

summarizeNames={'CorpusCallosum' 'unlabeled' 'OpticChi' 'ventricle' 'whiteMatter' 'pericallosal' 'insula'};

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

fprintf('\n categorical segmentation complete')
classificationOut
end