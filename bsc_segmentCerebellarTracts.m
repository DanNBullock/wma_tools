function [classificationOut] =bsc_segmentCerebellarTracts(wbfg, fsDir,experimentalBool,varargin)
% [classificationOut] =bsc_segmentCerebellarTracts(wbfg, fsDir,varargin)
%
% This function automatedly segments cerebellar
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.  Subsections may come later.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -experimentalBool: toggle for experimental tracts
% -varargin: priors from previous steps

% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University

%% parameter note & initialization

%create left/right lables
sideLabel={'left','right'};

%[categoryPrior] =bsc_streamlineCategoryPriors_v4(wbfg, fsDir,2)

categoryPrior=varargin{1};
%categoryPrior=categoryPrior{1};

%[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

cbROINums=[7 8;46 47];
thalamusLut=[10 49];

[inflatedAtlas] =bsc_inflateLabels(fsDir,2);

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

LeftCBRoi=bsc_roiFromAtlasNums(inflatedAtlas,cbROINums(1,:) ,1);
RightCBRoi=bsc_roiFromAtlasNums(inflatedAtlas,cbROINums(2,:) ,1);
[~, leftCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{LeftCBRoi}], {'endpoints'}, 'dud');
[~, rightCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{RightCBRoi}], {'endpoints'}, 'dud');
spinoCebBool=or(categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_spinal_interHemi')),categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_spinal')));
SpineTop= bsc_planeFromROI_v2(16,'superior',atlasPath);
[~, SpineTopBool]=wma_SegmentFascicleFromConnectome(wbfg, [{SpineTop}], {'not'}, 'dud');

classificationOut=bsc_concatClassificationCriteria(classificationOut,'leftSpinoCerebellar',leftCebBool,spinoCebBool,SpineTopBool);
classificationOut=bsc_concatClassificationCriteria(classificationOut,'rightSpinoCerebellar',rightCebBool,spinoCebBool,SpineTopBool);
SpineLimit= bsc_planeFromROI_v2(16,'anterior',atlasPath);
[~, posteriorStreams]=wma_SegmentFascicleFromConnectome(wbfg, [{SpineLimit}], {'not'}, 'dud');

for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% Ipsilateral connections
    CBRoi=bsc_roiFromAtlasNums(inflatedAtlas,cbROINums(leftright,:) ,1);
    ThalROI=bsc_roiFromAtlasNums(inflatedAtlas,thalamusLut(leftright) ,1);
    antCebSplit= bsc_planeFromROI_v2(thalamusLut(leftright),'anterior',atlasPath);
    
    %motorCerebellum
    frontoMotorLimit= bsc_planeFromROI_v2(170+sidenum, 'anterior',inflatedAtlas);
    frontMedROI=bsc_roiFromAtlasNums(inflatedAtlas,[116]+sidenum ,1);
    frontoMedMotorROI=bsc_modifyROI_v2(inflatedAtlas,frontMedROI, frontoMotorLimit, 'posterior');
    MotorROI=bsc_roiFromAtlasNums(inflatedAtlas,[168 128 146 129 170 103 ]+sidenum ,1);
    MotorROI=bsc_mergeROIs(frontoMedMotorROI,MotorROI);
    
    [~, AnterioFrontoBool]=wma_SegmentFascicleFromConnectome(wbfg, [{CBRoi} {antCebSplit}], {'endpoints','and'}, 'dud');
    [~, middleFrontoBool]=wma_SegmentFascicleFromConnectome(wbfg, [{CBRoi} {antCebSplit}], {'endpoints','not'}, 'dud');
    [~, thalCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{CBRoi} {ThalROI} {SpineLimit}], {'endpoints','endpoints','not'}, 'dud');
    [~, motorCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{CBRoi} {MotorROI}], {'endpoints','endpoints'}, 'dud');
    [~, thisCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{CBRoi}], {'endpoints'}, 'dud');
    
    motorCebBool=motorCebBool&or(categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_frontal')),categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_parietal')));
    AnterioFrontoBool=AnterioFrontoBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_frontal'));
    middleFrontoBool=~motorCebBool&middleFrontoBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_frontal'));
    thalCebBool=thalCebBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_subcortical'));
    occipitoCebBool=posteriorStreams&thisCebBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_occipital'));
    parietoCebBool=posteriorStreams&thisCebBool&~motorCebBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_parietal'));
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'MotorCerebellar'),motorCebBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},' AnterioFrontoCerebellar'),AnterioFrontoBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'MiddleFrontoBoolCerebellar'),middleFrontoBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ThalamicoCerebellar'),thalCebBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'OccipitoCerebellar'),occipitoCebBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ParietoCerebellar'),parietoCebBool);
    
    if experimentalBool
         classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'MiddleFrontoBoolCerebellar'),middleFrontoBool);
    end

    %% Contralateral connections
    if leftright==1
        otherCBRoi=bsc_roiFromAtlasNums(inflatedAtlas,cbROINums(2,:) ,1);
        otherWM=bsc_roiFromAtlasNums(atlasPath,41 ,1);
    else
        otherCBRoi=bsc_roiFromAtlasNums(inflatedAtlas,cbROINums(1,:) ,1);
         otherWM=bsc_roiFromAtlasNums(atlasPath,2 ,1);
    end
    
    [~,notThese]=wma_SegmentFascicleFromConnectome(wbfg, [{otherWM}], {'and'}, 'dud');
  
    
    [~, contraAnterioFrontoBool]=wma_SegmentFascicleFromConnectome(wbfg, [{otherCBRoi} {antCebSplit}], {'endpoints','and'}, 'dud');
    [~, contramiddleFrontoBool]=wma_SegmentFascicleFromConnectome(wbfg, [{otherCBRoi} {antCebSplit}], {'endpoints','not'}, 'dud');
    [~, contrathalCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{otherCBRoi} {ThalROI}, {SpineLimit}], {'endpoints','endpoints','not'}, 'dud');
    [~, contramotorCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{otherCBRoi} {MotorROI}], {'endpoints','endpoints'}, 'dud');
    [~, contrathisCebBool]=wma_SegmentFascicleFromConnectome(wbfg, [{otherCBRoi}], {'endpoints'}, 'dud');
    
    
    contramotorCebBool=~notThese&contramotorCebBool&or(categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_frontal_interHemi')),categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_parietal_interHemi')));
    contraAnterioFrontoBool=~notThese&contraAnterioFrontoBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_frontal_interHemi'));
    contramiddleFrontoBool=~notThese&~contramotorCebBool&contramiddleFrontoBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_frontal_interHemi'));
    contrathalCebBool=~notThese&contrathalCebBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_subcortical_interHemi'));
    contraoccipitoCebBool=~notThese&posteriorStreams&contrathisCebBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_occipital_interHemi'));
    contraparietoCebBool=~notThese&posteriorStreams&contrathisCebBool&~contramotorCebBool&categoryPrior.index==find(strcmp(categoryPrior.names,'cerebellum_to_parietal_interHemi'));
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ContraMotorCerebellar'),contramotorCebBool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ContraAnterioFrontoCerebellar'),contraAnterioFrontoBool);
     classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ContraOccipitoCerebellar'),contraoccipitoCebBool);
     classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ContraParietoCerebellar'),contraparietoCebBool);
     
     if experimentalBool
         classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ContraMiddleFrontoBoolCerebellar'),contramiddleFrontoBool);
         classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'ContraThalamicoCerebellar'),contrathalCebBool);
     end
     
end
end

