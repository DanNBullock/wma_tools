function [classificationOut] =bsc_segmentAntPostTracts_v2(wbfg, fsDir,varargin)
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

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

lentiLut=[12 13; 51 52];
palLut=[13;52];
thalLut=[10;49];
ventricleLut=[4;43];
wmLut=[2;41];
DCLut=[28;60];
hippLut=[17;53];
amigLut=[18;54];

subcort=[10 12 13 17 18; 49 51 52 53 54];

interHemiNot=bsc_makePlanarROI(atlasPath,0, 'x');

[classificationOut] =bsc_segmentCingulum_v3(wbfg, fsDir,categoryPrior);
cingulumBool=or(classificationOut.index==find(strcmp(classificationOut.names,'rightcingulum')),classificationOut.index==find(strcmp(classificationOut.names,'leftcingulum')));


%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    
    
    %%
    thalTop=bsc_planeFromROI_v2(thalLut(leftright), 'superior',atlasPath);
    thalPost=bsc_planeFromROI_v2(thalLut(leftright), 'posterior',atlasPath);
    amigPost=bsc_planeFromROI_v2(amigLut(leftright),'posterior',atlasPath);
    
    
    [~, UncSegBool]=wma_SegmentFascicleFromConnectome(wbfg, [{thalTop} {amigPost} {thalPost}], {'not','not','not'}, 'dud');
    
    frontoTemporalBool=or(bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_temporal')),bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_temporal_ufiber')));
    
    %categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),''),categoryPrior.names)))';
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Uncinate'),frontoTemporalBool,UncSegBool);
    
    %UNCINATE DONE ========================================================
    
    ccPostLimit=bsc_planeFromROI_v2(251, 'posterior',atlasPath);
    ccAntLimit=bsc_planeFromROI_v2(255, 'anterior',atlasPath);
    
    %carve out the area around and above the cc
    [postCCtopThal]=bsc_modifyROI_v2(atlasPath,ccPostLimit, thalTop, 'superior');
    [ccInterior1]=bsc_modifyROI_v2(atlasPath,thalTop, ccPostLimit, 'anterior');
    [ccInterior2]=bsc_modifyROI_v2(atlasPath,ccInterior1, ccAntLimit, 'posterior');
    [antCCtopThal]=bsc_modifyROI_v2(atlasPath,ccAntLimit, thalTop, 'superior');
    
    ccCarveOut=bsc_mergeROIs(postCCtopThal,ccInterior2);
    ccCarveOut=bsc_mergeROIs(ccCarveOut,antCCtopThal);
    
    antTempPlane=bsc_planeFromROI_v2(173+sidenum, 'anterior',atlasPath);
    infCCLimit=bsc_planeFromROI_v2(255, 'inferior',atlasPath);
    
    
    [infTempROI]=bsc_modifyROI_v2(atlasPath,antTempPlane, infCCLimit, 'inferior');
    
    [~, IFOFBool]=wma_SegmentFascicleFromConnectome(wbfg, [{infTempROI} {ccCarveOut}], {'and', 'not'}, 'dud');
    
    %[indexBool] = bsc_extractStreamIndByName(classification,tractName)
    frontoOccipitalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_occipital'));
    %frontoOccipitalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_occipital'),categoryPrior.names)))';
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'IFOF'),frontoOccipitalBool,IFOFBool);
    
    %IFOF DONE ========================================================
    %Create anatomical Rois
    latFisInf=bsc_planeFromROI_v2(141+sidenum, 'inferior',atlasPath);
    insPost=bsc_planeFromROI_v2(150+sidenum, 'posterior',atlasPath);
    tempTransVTop=bsc_planeFromROI_v2(133+sidenum, 'superior',atlasPath);
     
    postLatFisInf=bsc_modifyROI_v2(atlasPath,latFisInf, ccPostLimit, 'posterior');
       %neckArcAnd=bsc_modifyROI_v2(atlasPath,latFisInf, insPost, 'posterior');
    
    TopArcAnd=bsc_modifyROI_v2(atlasPath,insPost, tempTransVTop, 'superior');
    
    [~, arcBool]=wma_SegmentFascicleFromConnectome(wbfg, [{postLatFisInf} {TopArcAnd} ], {'and', 'and'}, 'dud');
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Arc'),frontoTemporalBool,arcBool,~cingulumBool);
    
    %Arcuate segmentation complete========================================
    
    %[indexBool] = bsc_extractStreamIndByName(classification,tractName)
    %parietoFrontalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_parietal'),categoryPrior.names)))';
    parietoFrontalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel(leftright),'frontal_to_parietal'));
    
    ccMidLimit=bsc_planeFromROI_v2(252, 'anterior',atlasPath);
        
    slf12exclude= bsc_modifyROI_v2(atlasPath,ccInterior2, ccMidLimit, 'posterior');

    %slf12ROI=bsc_roiFromAtlasNums(atlasPath,[115 114 116 154 155 153]+sidenum,1);
    
    [SFL3Intersection] = bsc_MultiIntersectROIs(atlasPath,19,112+sidenum, 150+sidenum );
    
        palAnt=bsc_planeFromROI_v2(palLut(leftright), 'anterior',atlasPath);
        frontSinfLimit=bsc_planeFromROI_v2(155+sidenum, 'inferior',atlasPath);
        
        slf3Fix= bsc_modifyROI_v2(atlasPath,palAnt, frontSinfLimit, 'superior');
    
    [~, SLF12Bool]=wma_SegmentFascicleFromConnectome(wbfg, [ {slf12exclude} {TopArcAnd}], {'not','and'}, 'dud');
    
    [~, SLF3Bool]=wma_SegmentFascicleFromConnectome(wbfg, [{SFL3Intersection} {postLatFisInf} {slf3Fix}], {'endpoints','not','not'}, 'dud');
    

    
    SLF12Bool=SLF12Bool&parietoFrontalBool&~cingulumBool&~IFOFBool;
    SLF3Bool=SLF3Bool&parietoFrontalBool&~IFOFBool;
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF1And2'),SLF12Bool);
    %classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF2'),SLF2Bool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF3'),SLF3Bool);
    
    
    
end


end