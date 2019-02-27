function [classificationOut] =bsc_streamlineCategoryPriors_v6(wbfg, fsDir,inflateITer)
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
leftROIS=[[101:1:175]+11000 26  17 18 7 8 10:13];
rightROIS=[[101:1:175]+12000 46 47 49:54 58];

subcorticalROIS=[10:13 17:20 26 58 27 49:56 59 ];
spineROIS=[16 28 60];
cerebellumROIS=[8 47 7 46 ];
ventricleROIS=[31 63 11 50 4 43 14 24 15 44 5 62 30 80 72 ];
wmROIS=[41 2];
ccROIS=[251:255];
unknownROIS=[0 2000 1000 77:82 24];
OpticCROI=[85];

FrontalROIs=[[124 148 118 165 101 154 105 115 154 155 115 170 129 146 153 ...
    164 106 116 108 131 171 112 150 104 169 114 113 116 107 163 139 132 140]+11000 [124 148 118 165 101 154 105 115 154 155 115 170 129 146 153 ...
    164 106 116 108 131 171 112 150 104 169 114 113 116 107 163 139 132 140]+12000] ;

TemporalROIs=[[144 134 138 137 173 174 135 175 121 151 123 162 133]+11000 [144 134 138 137 173 174 135 175 121 151 123 162 133]+12000];

OccipitalROI=[[120 119 111 158 166 143 145 159 152 122 162 161 121 160 102]+11000 [120 119 111 158 166 143 145 159 152 122 162 161 121 160 102]+12000];

ParietalROI=[[157 127 168 136 126 125 156 128 141 172 147 109 103 130 110]+11000 [157 127 168 136 126 125 156 128 141 172 147 109 103 130 110]+12000];

pericROI=[[167]+11000 [167]+12000];

insulaROI=[[117 149]+11000 [117 149]+12000];


endpoints1=zeros(3,length(allStreams));
endpoints2=zeros(3,length(allStreams));

for icategories=1:length(allStreams)
    curStream=allStreams{icategories};
    endpoints1(:,icategories)=curStream(:,1);
    endpoints2(:,icategories)=curStream(:,end);
end

if inflateITer>0
    [inflatedAtlas] =bsc_inflateLabels(fsDir,inflateITer);
else
    inflatedAtlas=niftiRead(atlasPath);
end

[endpoints1Identity] =bsc_atlasROINumsFromCoords_v3(inflatedAtlas,endpoints1,'acpc');
[endpoints2Identity] =bsc_atlasROINumsFromCoords_v3(inflatedAtlas,endpoints2,'acpc');
frpintf('\n endpoint identities determined')

excludeBool=zeros(1,length(allStreams));
includeBool=excludeBool;
LeftBool=excludeBool;
RightBool=excludeBool;
interHemiBool=excludeBool;
validUIndexes=excludeBool;
singleLeftBoolproto=excludeBool;
singleRightBoolproto=excludeBool;
interhemiFlag=excludeBool;
termination2=cell(1,length(allStreams));
termination1=cell(1,length(allStreams));
streamName=termination1;

[superficialClassification] =bsc_segmentSuperficialFibers(wbfg, fsDir);
frpintf('\n superficial fibers identified')

validSideROI= [leftROIS rightROIS] ;
excludeSideROI=[unknownROIS pericROI ccROIS OpticCROI wmROIS spineROIS ventricleROIS];
  
for iStreams=1:length(allStreams)
% disagreeBool(iStreams)=or(and(LeftBool(iStreams),and(ismember(endpoints2Identity(iStreams),leftROIS),ismember(endpoints1Identity(iStreams),leftROIS))),and(RightBool(iStreams),and(ismember(endpoints2Identity(iStreams),rightROIS),ismember(endpoints1Identity(iStreams),rightROIS))))  ;      
excludeBool(iStreams)=or(ismember(endpoints2Identity(iStreams),excludeSideROI),ismember(endpoints1Identity(iStreams),excludeSideROI));
includeBool(iStreams)=or(ismember(endpoints2Identity(iStreams),validSideROI),ismember(endpoints1Identity(iStreams),validSideROI));
validUIndexes(iStreams)=or(ismember(endpoints2Identity(iStreams),greyMatterROIS),ismember(endpoints1Identity(iStreams),greyMatterROIS));
LeftBool(iStreams)=and(ismember(endpoints2Identity(iStreams),leftROIS),ismember(endpoints1Identity(iStreams),leftROIS));
RightBool(iStreams)=and(ismember(endpoints2Identity(iStreams),rightROIS),ismember(endpoints1Identity(iStreams),rightROIS));
interHemiBool(iStreams)=or(and(ismember(endpoints2Identity(iStreams),leftROIS),ismember(endpoints1Identity(iStreams),rightROIS)),and(ismember(endpoints2Identity(iStreams),rightROIS),ismember(endpoints1Identity(iStreams),leftROIS)));

singleLeftBoolproto(iStreams)=xor(ismember(endpoints2Identity(iStreams),leftROIS),ismember(endpoints1Identity(iStreams),leftROIS));
singleRightBoolproto(iStreams)=xor(ismember(endpoints2Identity(iStreams),rightROIS),ismember(endpoints1Identity(iStreams),rightROIS));




if     ~isempty(find(endpoints1Identity(iStreams)==FrontalROIs, 1))
    termination1{iStreams}='frontal';
elseif ~isempty(find(endpoints1Identity(iStreams)==TemporalROIs, 1))
    termination1{iStreams}='temporal';
elseif ~isempty(find(endpoints1Identity(iStreams)==OccipitalROI, 1))
    termination1{iStreams}='occipital';
elseif ~isempty(find(endpoints1Identity(iStreams)==ParietalROI, 1))
    termination1{iStreams}='parietal';
elseif ~isempty(find(endpoints1Identity(iStreams)==subcorticalROIS, 1))
    termination1{iStreams}='subcortical';
elseif ~isempty(find(endpoints1Identity(iStreams)==spineROIS, 1))
    termination1{iStreams}='spinal';
elseif ~isempty(find(endpoints1Identity(iStreams)==insulaROI, 1))
    termination1{iStreams}='insula';
elseif ~isempty(find(endpoints1Identity(iStreams)==cerebellumROIS, 1))
    termination1{iStreams}='cerebellum';
elseif ~isempty(find(endpoints1Identity(iStreams)==ccROIS, 1))
    termination1{iStreams}='CorpusCallosum';
    %false positives
elseif ~isempty(find(endpoints1Identity(iStreams)==ventricleROIS, 1))
    termination1{iStreams}='ventricle';
elseif ~isempty(find(endpoints1Identity(iStreams)==unknownROIS, 1))
    termination1{iStreams}='unlabeled';
elseif ~isempty(find(endpoints1Identity(iStreams)==wmROIS, 1))
    termination1{iStreams}='whiteMatter';
elseif ~isempty(find(endpoints1Identity(iStreams)==pericROI, 1))
    termination1{iStreams}='pericallosal';
elseif ~isempty(find(endpoints1Identity(iStreams)==OpticCROI, 1))
    termination1{iStreams}='OpticChi';
end

if     ~isempty(find(endpoints2Identity(iStreams)==FrontalROIs, 1))
    termination2{iStreams}='frontal';
elseif ~isempty(find(endpoints2Identity(iStreams)==TemporalROIs, 1))
    termination2{iStreams}='temporal';
elseif ~isempty(find(endpoints2Identity(iStreams)==OccipitalROI, 1))
    termination2{iStreams}='occipital';
elseif ~isempty(find(endpoints2Identity(iStreams)==ParietalROI, 1))
    termination2{iStreams}='parietal';
elseif ~isempty(find(endpoints2Identity(iStreams)==subcorticalROIS, 1))
    termination2{iStreams}='subcortical';
elseif ~isempty(find(endpoints2Identity(iStreams)==spineROIS, 1))
    termination2{iStreams}='spinal';
elseif ~isempty(find(endpoints2Identity(iStreams)==insulaROI, 1))
    termination2{iStreams}='insula';
elseif ~isempty(find(endpoints2Identity(iStreams)==cerebellumROIS, 1))
    termination2{iStreams}='cerebellum';
    %false positives
elseif ~isempty(find(endpoints2Identity(iStreams)==pericROI, 1))
    termination2{iStreams}='pericallosal';
elseif ~isempty(find(endpoints2Identity(iStreams)==ventricleROIS, 1))
    termination2{iStreams}='ventricle';
elseif ~isempty(find(endpoints2Identity(iStreams)==unknownROIS, 1))
    termination2{iStreams}='unlabeled';
elseif ~isempty(find(endpoints2Identity(iStreams)==wmROIS, 1))
    termination2{iStreams}='whiteMatter';
elseif ~isempty(find(endpoints2Identity(iStreams)==ccROIS, 1))
    termination2{iStreams}='CorpusCallosum';
elseif ~isempty(find(endpoints2Identity(iStreams)==OpticCROI, 1))
    termination2{iStreams}='OpticChi';
end


if ~or(isempty(termination1{iStreams}),isempty(termination2{iStreams}))
    terminationNames=sort({termination1{iStreams} termination2{iStreams}});
else
    endpoints1Identity(iStreams)
    endpoints2Identity(iStreams)
    error('streamline identity unaccounted for')
end

fprintf('\n endpoint categories determined')

%hierarchy of categories here
interhemiFlag(iStreams)=interHemiBool(iStreams)&includeBool(iStreams);
if interhemiFlag(iStreams)
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2},'_interHemi');
else
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2});
end
if superficialClassification.index(iStreams)>0&validUIndexes(iStreams)
    streamName{iStreams}=strcat(terminationNames{1},'_to_',terminationNames{2},'_ufiber');
end

if or(LeftBool(iStreams),singleLeftBoolproto(iStreams))&includeBool(iStreams)&~interhemiFlag(iStreams)
    streamName{iStreams}=strcat('left',streamName{iStreams});
elseif or(RightBool(iStreams),singleRightBoolproto(iStreams))&includeBool(iStreams)&~interhemiFlag(iStreams)
    streamName{iStreams}=strcat('right',streamName{iStreams});
end

end

uniqueNames=unique(streamName);



summarizeNames={'CorpusCallosum' 'unlabeled' 'OpticChi' 'ventricle' 'whiteMatter' 'pericallosal'};

for icategories=1:length(uniqueNames)
    icategories
    if contains(uniqueNames{icategories},summarizeNames)
        for isummary=1:length(summarizeNames)
            summaryIndex(isummary)=contains(uniqueNames{icategories},summarizeNames{isummary});
        end
        summaryIndexSingle=find(summaryIndex);
        classificationOut=bsc_concatClassificationCriteria(classificationOut,summarizeNames{summaryIndexSingle(1)},contains(streamName,uniqueNames{icategories}));
        clear summaryIndex
    else
        classificationOut=bsc_concatClassificationCriteria(classificationOut,uniqueNames{icategories},ismember(streamName,uniqueNames{icategories}));
    end
end

classificationOut = wma_resortClassificationStruc(classificationOut);

fprintf('\n categorical segmentation complete')
end