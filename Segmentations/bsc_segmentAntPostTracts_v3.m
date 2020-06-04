function [classificationOut] =bsc_segmentAntPostTracts_v3(wbfg,atlas,varargin)
% [classificationOut] =bsc_segmentArc_Cingulum(wbfg, fsDir)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -varargin: priors from previous steps

% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University

%% parameter note & initialization

%create left/right lables
sideLabel={'left','right'};

categoryPrior=varargin{1};
effPrior=varargin{2};

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

%atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

lentiLut=[12 13; 51 52];
palLut=[13;52];
thalLut=[10;49];
ventricleLut=[4;43];
wmLut=[2;41];
DCLut=[28;60];
hippLut=[17;53];
amigLut=[18;54];

subcort=[10 12 13 17 18; 49 51 52 53 54];

interHemiNot=bsc_makePlanarROI(atlas,0, 'x');

%do it here as a prior
[classificationOut] =bsc_segmentCingulum_v3(wbfg,atlas,categoryPrior);
cingulumBool=or(classificationOut.index==find(strcmp(classificationOut.names,'rightcingulum')),classificationOut.index==find(strcmp(classificationOut.names,'leftcingulum')));

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    thalTop=bsc_planeFromROI_v2(thalLut(leftright), 'superior',atlas);
    thalPost=bsc_planeFromROI_v2(thalLut(leftright), 'posterior',atlas);
    amigPost=bsc_planeFromROI_v2(amigLut(leftright),'posterior',atlas);
    
    %ensure endpoints are anterior of specified planes
    uncEndpointCriteria=bsc_applyEndpointCriteria(wbfg,thalTop,'inferior','both',amigPost,'anterior','both');
    
    [~, UncSegBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalTop} {amigPost} {thalPost}], {'not','not','not'}, 'dud');
    
    frontoTemporalBool=or(bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_temporal')), bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_temporal_ufiber')));    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Uncinate'),frontoTemporalBool,UncSegBool,uncEndpointCriteria);
    
    %UNCINATE DONE ========================================================
    
    ccPostLimit=bsc_planeFromROI_v2(251, 'posterior',atlas);
    ccAntLimit=bsc_planeFromROI_v2(255, 'anterior',atlas);
    
    %carve out the area around and above the cc
    [postCCtopThal]=bsc_modifyROI_v2(atlas,ccPostLimit, thalTop, 'superior');
    [ccInterior1]=bsc_modifyROI_v2(atlas,thalTop, ccPostLimit, 'anterior');
    [ccInterior2]=bsc_modifyROI_v2(atlas,ccInterior1, ccAntLimit, 'posterior');
    [antCCtopThal]=bsc_modifyROI_v2(atlas,ccAntLimit, thalTop, 'superior');
    
    ccCarveOut=bsc_mergeROIs(postCCtopThal,ccInterior2);
    ccCarveOut=bsc_mergeROIs(ccCarveOut,antCCtopThal);
    
    antTempPlane=bsc_planeFromROI_v2(173+sidenum, 'anterior',atlas);
    infCCLimit=bsc_planeFromROI_v2(255, 'inferior',atlas);
    
    
    [infTempROI]=bsc_modifyROI_v2(atlas,antTempPlane, infCCLimit, 'inferior');
    
    [~, IFOFBool]=wma_SegmentFascicleFromConnectome(wbfg, [{infTempROI} {ccCarveOut}], {'and', 'not'}, 'dud');
    
    %[indexBool] = bsc_extractStreamIndByName(classification,tractName)
    frontoOccipitalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_occipital'));
    %frontoOccipitalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_occipital'),categoryPrior.names)))';
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'IFOF'),frontoOccipitalBool,IFOFBool);
    
    %IFOF DONE ========================================================
    %Create anatomical Rois
    
    postCingNot=bsc_roiFromAtlasNums(atlas,[108]+sidenum ,5);
    latFisInf=bsc_planeFromROI_v2(141+sidenum, 'inferior',atlas);
    insPost=bsc_planeFromROI_v2(150+sidenum, 'posterior',atlas);
    tempTransVTop=bsc_planeFromROI_v2(133+sidenum, 'superior',atlas);
     
    postLatFisInf=bsc_modifyROI_v2(atlas,latFisInf, ccPostLimit, 'posterior');

    TopArcAnd=bsc_modifyROI_v2(atlas,insPost, tempTransVTop, 'superior');
    
    [~, arcBool]=wma_SegmentFascicleFromConnectome(wbfg, [{postLatFisInf} {TopArcAnd} {postCingNot} ], {'and', 'and','not'}, 'dud');
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Arc'),frontoTemporalBool,arcBool,~cingulumBool);
    
    %Arcuate segmentation complete========================================
    
    %[indexBool] = bsc_extractStreamIndByName(classification,tractName)
    %parietoFrontalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_parietal'),categoryPrior.names)))';
    parietoFrontalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_parietal'));
    
    ccMidLimit=bsc_planeFromROI_v2(252, 'anterior',atlas);
        
    slf12exclude= bsc_modifyROI_v2(atlas,ccInterior2, ccMidLimit, 'posterior');
    
    [interiorWallBool]=bsc_endpointAtlasCriteria(wbfg,atlas,[116 109 108  107 167 147 172 130 110 106]+sidenum,'either');
   
    [SFL3Intersection] = bsc_MultiIntersectROIs(atlas,19,112+sidenum, 150+sidenum );
    
    palAnt=bsc_planeFromROI_v2(palLut(leftright), 'anterior',atlas);
    frontSinfLimit=bsc_planeFromROI_v2(155+sidenum, 'inferior',atlas);
    
    slf3Fix= bsc_modifyROI_v2(atlas,palAnt, frontSinfLimit, 'superior');
    
    [~, SLF12segBool]=wma_SegmentFascicleFromConnectome(wbfg, [ {slf12exclude} {TopArcAnd} ], {'not','and'}, 'dud');
    
    [~, SLF3Bool]=wma_SegmentFascicleFromConnectome(wbfg, [{SFL3Intersection} {postLatFisInf} {slf3Fix} ], {'endpoints','not','not'}, 'dud');
    
    
    SLF12Bool=SLF12segBool&parietoFrontalBool&~cingulumBool&~IFOFBool&~interiorWallBool;
    SLF3Bool=SLF3Bool&parietoFrontalBool&~IFOFBool&~interiorWallBool;
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF1And2'),SLF12Bool);
    %classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF2'),SLF2Bool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF3'),SLF3Bool);
    
    
    %%
    
    
end


end