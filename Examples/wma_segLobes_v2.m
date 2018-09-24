function [classification] =wma_segLobes_v2(feORwbfg, fsDir)
% [classification] =wma_segLobes(wbfg, fsDir)
% This function automatedly segments the intralobal and interlobal fibers
% of the brain;
%
% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% classification:  a classification structure with .name and .indexes
% fields
%
% (C) Daniel Bullock, 2018, Indiana University
%
%% parameter note & initialization

% extract fg if necessary
[wbfg, ~] = bsc_LoadAndParseFiberStructure(feORwbfg);

% find midpoints
for iFibers=1:length(wbfg.fibers)
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
end

%set Left/Right parameters
LeftStreams=midpoints(:,1)<0;
sideLabel={'left','right'};

%initialize classification structure
classification=[];
classification.index=zeros(length(wbfg.fibers),1);

%set cerebellar ROIs (these do not follow the same numerical patterns as
%cortical ROIs)
cerebellumROI=[8 7; 46 47];

%% segmentation
for leftright= [1,2]
    %Set numerical ROI hemisphere indicator (i.e. 12000 = right, 11000 =
    %left)
    sidenum=10000+leftright*1000;
    
    %Create lobe ROIs;
    FrontalROI=bsc_roiFromFSnums(fsDir,[124 148 118 165 101 154 105 115 154 155 170 129 146 153 ...
        164 106 116 108 131 171 112 150 104 169 114 113 103 107 117 132 139 140 142 163]+sidenum ,0);
    TemporalROI=bsc_roiFromFSnums(fsDir,[144 134 138 137 173 174 135 175 121 151 123 162 133 149 ]+sidenum ,0);
    OccipitalROI=bsc_roiFromFSnums(fsDir,[120 119 111 158 166 143 145 159 152 122 162 161 102 160]+sidenum,0);
    ParietalROI=bsc_roiFromFSnums(fsDir,[157 127 168 136 126 125 156 128 141 172 147 109 110 130]+sidenum,0);
    CerebellarROI=bsc_roiFromFSnums(fsDir,cerebellumROI(leftright,:),0);
    
    %% WMA toolkit demonstration:
    %Note, we did not include the pericollosal ROI (11167/12167) due to
    %it's traversing multiple lobes.  Here is an example of how to deal with
    %that:
    
    %looking at freesurfer's interactive visualization we can note that the
    %anterior of the subparietal sulcus ROI (11172/12172) can serve as a good
    %demarcation of where the frontal portion of the pericollosal ROI would
    %end.  Lets go ahead and extract a planar roi from it.
    [antSubPar] = bsc_planeFromROI(172+sidenum, 'anterior',fsDir);
    %this will serve as the cutoff point for the pericollosal ROI, we want
    %everything anterior of it to merge into the frontal ROI.
    [antPeriCall]=bsc_modifyROI(fsDir,167+sidenum, antSubPar, 'anterior');
    %we can then merge this ROI into the preexisting Frontal ROI
    [FrontalROI] = bsc_mergeROIs(FrontalROI, antPeriCall);
    %however, we would also like to add the corresponding parts of the
    %posterior pericollosal ROI to the respictive lobes.  Lets begin by
    %obtaining that portion of the ROI.
    [postPeriCall]=bsc_modifyROI(fsDir,167+sidenum, antSubPar, 'posterior');
    %this posterior portion must now be separated into parietal and
    %temporal components.  Again looking at the interactive freesurfer
    %visualization we note that the post-dorsal cingulum ROI (11109/12109)
    %demarkates the separation between the temporal and parietal components
    %of the previously made postPeriCall ROI.  We can thus sub-separate
    %this ROI at this point.  Lets begin by extracting the relevant planar
    %ROI
    [infPstDorCing] = bsc_planeFromROI(109+sidenum, 'inferior',fsDir);
    %using this ROI we can separate the postPeriCall ROI
    [superiorPostPeriCall]=bsc_modifyROI(fsDir,postPeriCall, infPstDorCing, 'superior');
    [inferiorPostPeriCall]=bsc_modifyROI(fsDir,postPeriCall, infPstDorCing, 'inferior');
    %merge these into the respective ROIs
    [TemporalROI] = bsc_mergeROIs(TemporalROI, inferiorPostPeriCall);
    [ParietalROI] = bsc_mergeROIs(ParietalROI, superiorPostPeriCall);
    %Now we can proceed with the rest of the segmentation
    
    %% continuing with the segmentation
    %Segment within lobe tracts
    [~, FrontalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {FrontalROI}], { 'both_endpoints' }, 'Frontal');
    [~, TemporalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {TemporalROI}], { 'both_endpoints' }, 'Temporal');
    [~, OccipitalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {OccipitalROI}], { 'both_endpoints' }, 'Occipital');
    [~, ParietalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {ParietalROI}], { 'both_endpoints' }, 'Parietal');
    [~, CerebellarIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {CerebellarROI}], { 'both_endpoints' }, 'Cerebellar');
    
    %Segment between lobe tracts
    [~, FrontalTemporalInd]=bsc_tractByEndpointROIs(wbfg,[{FrontalROI},{TemporalROI}]);
    [~, OccipitalTemporalInd]=bsc_tractByEndpointROIs(wbfg,[{OccipitalROI},{TemporalROI}]);
    [~, ParietalTemporalInd]=bsc_tractByEndpointROIs(wbfg,[{ParietalROI},{TemporalROI}]);
    [~, OccipitalParietalInd]=bsc_tractByEndpointROIs(wbfg,[{OccipitalROI},{ParietalROI}]);
    [~, OccipitalFrontalInd]=bsc_tractByEndpointROIs(wbfg,[{OccipitalROI},{FrontalROI}]);
    [~, FrontalParietalInd]=bsc_tractByEndpointROIs(wbfg,[{FrontalROI},{ParietalROI}]);
    
    %Assign streamlines to output classification structure
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'IntraFrontal'),FrontalIND);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'IntraTemporal'),TemporalIND);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'IntraOccipital'),OccipitalIND);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'IntraParietal'),ParietalIND);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'IntraCerebellar'),CerebellarIND);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'FrontoTemporal'),FrontalTemporalInd);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'OccipitoTemporal'),OccipitalTemporalInd);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'ParietoTemporal'),ParietalTemporalInd);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'OccipioParietal'),OccipitalParietalInd);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'FrontoOccipitallInd'),OccipitalFrontalInd);
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'FrontoParietal'),FrontalParietalInd);

end
end