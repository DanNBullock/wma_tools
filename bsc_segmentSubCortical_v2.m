function [classificationOut] =bsc_segmentSubCortical_v2(wbfg, fsDir,varargin)
%  [classificationOut] =bsc_segmentSubCortical(wbfg, fsDir,varargin)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -varargin{1}: categoryPrior, from bsc_streamlineCategoryPriors_v3 
% -varargin{2}: effPrior, from ConnectomeTestQ_v2
%
% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University
%% Begin code
%create left/right labels



categoryPrior=varargin{1};
effPrior=varargin{2};

atlasPath=fullfile(fsDir,'mri/aparc.a2009s+aseg.nii.gz');

sideLabel={'left','right'};

classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

lentiLut=[12 13; 51 52];
putLut=[12 ; 51 ];
palLut=[ 13; 52];
thalamusLut=[10 49];

[inflatedAtlas] =bsc_inflateLabels(fsDir,2);

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    FrontalROI=bsc_roiFromAtlasNums(inflatedAtlas,[124 148 118 165 101 154 105 115 154 155 170 129 146 153 ...
        164 106 116 108 131 171 112 150 104 169 114 113 103 107 117 132 139 140 142 163]+sidenum ,1);
    TemporalROI=bsc_roiFromAtlasNums(inflatedAtlas,[144 134 138 137 173 174 135 175 121 151 123 162 133 149 ]+sidenum ,1);
    OccipitalROI=bsc_roiFromAtlasNums(inflatedAtlas,[120 119 111 158 166 143 145 159 152 122 162 161 102 160]+sidenum,1);
    ParietalROI=bsc_roiFromAtlasNums(inflatedAtlas,[157 127 168 136 126 125 156 128 141 172 147 109 110 130]+sidenum,1);
    lentiROI=bsc_roiFromAtlasNums(inflatedAtlas,lentiLut(leftright,:),1);
    thalamicROI=bsc_roiFromAtlasNums(inflatedAtlas,thalamusLut(leftright),1);
    palROI=bsc_roiFromAtlasNums(inflatedAtlas,palLut(leftright),1);
    putROI=bsc_roiFromAtlasNums(inflatedAtlas,putLut(leftright),1);
    
    
    frontoMotorLimit= bsc_planeFromROI_v2(170+sidenum, 'anterior',inflatedAtlas);
    frontMedROI=bsc_roiFromAtlasNums(inflatedAtlas,[116]+sidenum ,1);
    frontoMedMotorROI=bsc_modifyROI_v2(inflatedAtlas,frontMedROI, frontoMotorLimit, 'posterior');
    MotorROI=bsc_roiFromAtlasNums(inflatedAtlas,[168 128 146 129 170 103 ]+sidenum ,1);
    
    
    MotorROI=bsc_mergeROIs(frontoMedMotorROI,MotorROI);
    %108?
    inferiorThalLimit= bsc_planeFromROI_v2(thalamicROI, 'inferior',atlasPath);
    anteriorThalLimit= bsc_planeFromROI_v2(thalamicROI, 'anterior',atlasPath);
    posteriorThalLimit= bsc_planeFromROI_v2(thalamicROI, 'posterior',atlasPath);
    superiorThalLimit= bsc_planeFromROI_v2(thalamicROI, 'superior',atlasPath);
    lateralThalamicLimit=bsc_planeFromROI_v2(thalamicROI, 'lateral',atlasPath);
    posteriorLentiLimit=bsc_planeFromROI_v2(lentiLut(leftright,:), 'posterior',atlasPath);
    antPalLimit=bsc_planeFromROI_v2(palLut(leftright,:), 'anterior',atlasPath);
    postPalLimit=bsc_planeFromROI_v2(palLut(leftright,:), 'posterior',atlasPath);
    
    antSpineLimit=bsc_planeFromROI_v2(16, 'anterior',atlasPath);
    postPutExclude=bsc_modifyROI_v2(inflatedAtlas,putROI, antSpineLimit, 'posterior');
    
    
    superiorFrontalLimit=bsc_planeFromROI_v2(112+sidenum, 'superior',atlasPath);
    anteriorFrontalLimit=bsc_planeFromROI_v2(155+sidenum, 'anterior',atlasPath);
    anteriorSuperiorFrontalLimit=bsc_modifyROI_v2(inflatedAtlas,superiorFrontalLimit, anteriorFrontalLimit, 'posterior');
    
    postPeriC= bsc_planeFromROI_v2(167+sidenum, 'posterior',atlasPath);
    
    posteriorInferiorExclude=bsc_modifyROI_v2(inflatedAtlas,inferiorThalLimit, anteriorThalLimit, 'posterior');
    posteriorSuperiorExclude=bsc_modifyROI_v2(inflatedAtlas,superiorThalLimit, anteriorThalLimit, 'posterior');
    posteriorLateralExclude=bsc_modifyROI_v2(inflatedAtlas,lateralThalamicLimit, anteriorThalLimit, 'posterior');
    anteriorLateralExclude=bsc_modifyROI_v2(inflatedAtlas, antPalLimit,lateralThalamicLimit, 'lateral');
    middleLateralExclude=bsc_modifyROI_v2(inflatedAtlas,posteriorLentiLimit, lateralThalamicLimit, 'medial');
    postPalLateralExclude=bsc_modifyROI_v2(inflatedAtlas, lateralThalamicLimit,postPalLimit, 'posterior');
    
    
    anteriorMedialExclude=bsc_modifyROI_v2(inflatedAtlas,anteriorThalLimit,lateralThalamicLimit, 'medial');
    %anteriorMedialExclude=bsc_modifyROI_v2(inflatedAtlas,anteriorThalLimit,lateralThalamicLimit, 'medial');
    
    [~, FrontoMotorBool]=wma_SegmentFascicleFromConnectome(wbfg, [{MotorROI}], {'endpoints'}, 'dud');
    [~, frontalBool]=wma_SegmentFascicleFromConnectome(wbfg, [{FrontalROI} {posteriorInferiorExclude} {posteriorSuperiorExclude} {posteriorLateralExclude}], {'endpoints','not','not','not'}, 'dud');
    [~, temporalBool]=wma_SegmentFascicleFromConnectome(wbfg, [{TemporalROI}], {'endpoints'}, 'dud');
    [~, occipitalBool]=wma_SegmentFascicleFromConnectome(wbfg, [{OccipitalROI}], {'endpoints'}, 'dud');
    [~, parietalBool]=wma_SegmentFascicleFromConnectome(wbfg, [{ParietalROI} {posteriorInferiorExclude} {anteriorThalLimit}], {'endpoints', 'not', 'not'}, 'dud');
    %create different versions of this
    [~, frontoThalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} {posteriorLateralExclude} {posteriorInferiorExclude} {middleLateralExclude} {anteriorSuperiorFrontalLimit}], {'endpoints', 'not', 'not' ,'not','not'}, 'dud');
    [~, temporoThalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} {postPeriC} {anteriorLateralExclude}], {'endpoints','not','not'}, 'dud');
    [~, parietoThalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} {antPalLimit} {inferiorThalLimit}], {'endpoints','not', 'not'}, 'dud');
    [~, lentiThalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} {palROI} {postPutExclude}], {'endpoints','endpoints','not'}, 'dud');
    [~, motorTthalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} {inferiorThalLimit} {anteriorLateralExclude}], {'endpoints','not','not'}, 'dud');
    [~, spinoThalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} {superiorThalLimit} {anteriorThalLimit}], {'endpoints','not','not'}, 'dud');
    
    [~, thalamicBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI} ], {'endpoints'}, 'dud');
    [~, lentiBool]=wma_SegmentFascicleFromConnectome(wbfg, [{lentiROI} {thalamicROI}], {'endpoints', 'endpoints'}, 'dud');
    
    
    % MAKE FRONTAL ENDPOINT CRITERIA
    frontoEndpointBool=bsc_applyEndpointCriteria(wbfg,anteriorThalLimit,'anterior','one');
   
    frontalSubcorticalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_subcortical'),categoryPrior.names)))';
    parietalSubcorticalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'parietal_to_subcortical'),categoryPrior.names)))';
    temporalSubcorticalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'subcortical_to_temporal'),categoryPrior.names)))';
    spinalSubcorticalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'spinal_to_subcortical'),categoryPrior.names)))';
     intraSubcorticalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'subcortical_to_subcortical'),categoryPrior.names)))';
    
    motorBoolPrior=or(parietalSubcorticalBool,frontalSubcorticalBool);
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'frontoThalamic'),frontoThalamicBool,frontoEndpointBool,frontalBool,~FrontoMotorBool,frontalSubcorticalBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'temporoThalamic'),temporoThalamicBool,temporalBool,temporalSubcorticalBool);
    % no need to duplicate the
    %classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'occipitalThalamic'),thalamicBool,occipitalBool,categoryPrior.index==find(strcmp(categoryPrior.names,'cortex_to_subcortical')));
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'parietoThalamic'),parietoThalamicBool,parietalBool,parietalSubcorticalBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'lentiThalamic'),lentiThalamicBool,lentiBool,intraSubcorticalBool);%
     classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'motorThalamic'),motorTthalamicBool,~frontalBool,FrontoMotorBool,motorBoolPrior);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'spinoThalamic'),spinoThalamicBool,spinalSubcorticalBool);
    
    
end
