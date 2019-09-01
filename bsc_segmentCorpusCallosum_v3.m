function [classificationOut] =bsc_segmentCorpusCallosum_v3(wbfg,atlas,experimentalBool,varargin)
%
%[RightILF, RightILFIndexes, LeftILF, LeftILFIndexes, RightMdLFspl, RightMdLFsplIndexes, LeftMdLFspl, LeftMdLFsplIndexes,...
%    RightMdLFang, RightMdLFangIndexes, LeftMdLFang, LeftMdLFangIndexes] =bsc_segmentMdLF_ILF(wbfg, fsDir)
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

categoryPrior=varargin{1};

%atlas=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

fullCCROI=bsc_roiFromAtlasNums(atlas,[255 254 253 252 251],13);
[~, excludeCCBool]=wma_SegmentFascicleFromConnectome(wbfg, [{fullCCROI}], {'endpoints'}, 'dud');
ventricleLut=[4;43];
ventricleROI=bsc_roiFromAtlasNums(atlas,[ventricleLut]',7);
[~, excludeVentriBool]=wma_SegmentFascicleFromConnectome(wbfg, [{ventricleROI}], {'endpoints'}, 'dud');




%create left/right lables
sideLabel={'left','right'};

classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);


lentiLut=[12 13; 51 52];
palLut=[13;52];
thalLut=[10;49];

wmLut=[2;41];
DCLut=[28;60];

subcort=[10 12 13 17 18; 49 51 52 53 54];

%sidenum is basically a way of switching  between the left and right
%hemispheres of the brain in accordance with freesurfer's ROI
%numbering scheme. left = 1, right = 2

anteriorCCROI=bsc_roiFromAtlasNums(atlas,255,0);
CCBottom= bsc_planeFromROI_v2(255, 'inferior',atlas);
[antMidCCPlaneAnt] = bsc_planeFromROI_v2(254, 'anterior',atlas);
[ccExtremAntPlaneAnt] = bsc_planeFromROI_v2(255, 'anterior',atlas);
[periCROI] = bsc_planeFromROI_v2(167+[11000 12000], 'anterior',atlas);
frontMiddlePlane=bsc_planeFromROI_v2(106+[11000 12000], 'anterior',atlas);
frontMiddleTopPlane=bsc_planeFromROI_v2(154+[11000 12000], 'superior',atlas);
olfTop=bsc_planeFromROI_v2(164+[11000 12000], 'superior',atlas);
ccAntPlane=bsc_planeFromROI_v2(255, 'anterior',atlas);



forcepsMiniorBool=bsc_applyEndpointCriteria(wbfg, periCROI,'anterior','both',olfTop,'superior','both',frontMiddleTopPlane,'inferior','both');
%anteriorCCPlaneBool=bsc_applyEndpointCriteria(wbfg, ccAntPlane,'anterior','both');

[anteriorCC, anteriorCCBool]=wma_SegmentFascicleFromConnectome(wbfg, [{anteriorCCROI} {antMidCCPlaneAnt} {periCROI}], {'and','not','and'}, 'dud');

frontalCCBool=bsc_extractStreamIndByName(categoryPrior,'frontal_to_frontal_interHemi');
anteriorCC.fibers=wbfg.fibers(forcepsMiniorBool&anteriorCCBool&~excludeCCBool&~excludeVentriBool&frontalCCBool);

forcepsMinorBool=forcepsMiniorBool&anteriorCCBool&~excludeCCBool&~excludeVentriBool&frontalCCBool;
classificationOut=bsc_concatClassificationCriteria(classificationOut,'forcepsMinor',forcepsMinorBool);




  LatTempPostLimit=bsc_planeFromROI_v2(161+[11000 12000], 'posterior',atlas);
  posteriorCCInfLimit=bsc_planeFromROI_v2(251, 'inferior',atlas);
 postCCinferiorROIcut=bsc_modifyROI_v2(atlas,posteriorCCInfLimit, LatTempPostLimit, 'anterior');
 [ccPostMidpoints]=bsc_applyMidpointCriteria(wbfg, LatTempPostLimit,'anterior');


 [posteriorCC, posteriorCCBool]=wma_SegmentFascicleFromConnectome(wbfg, [{postCCinferiorROIcut}], {'not'}, 'dud');
 
 posteriorCCBool=or(bsc_extractStreamIndByName(categoryPrior,'occipital_to_occipital_interHemi'),bsc_extractStreamIndByName(categoryPrior,'MaskFailure'));

forcepsMajorBool=ccPostMidpoints&posteriorCCBool&posteriorCCBool;
classificationOut=bsc_concatClassificationCriteria(classificationOut,'forcepsMajor',forcepsMajorBool);




parietalCCposteriorLimit=bsc_planeFromROI_v2(253, 'posterior',atlas);
parietalCCinferiorLimit=bsc_planeFromROI_v2(251, 'inferior',atlas);
parietalCCminTopLimit=bsc_planeFromROI_v2(252, 'superior',atlas);
latPutLeft=bsc_planeFromROI_v2(12, 'lateral',atlas);
latPutRight=bsc_planeFromROI_v2(51, 'lateral',atlas);
parietalEndpointsBool=bsc_applyEndpointCriteria(wbfg, parietalCCminTopLimit,'superior','both');

[parietalCC, parietalCCBool]=wma_SegmentFascicleFromConnectome(wbfg, [{parietalCCposteriorLimit parietalCCinferiorLimit} {latPutLeft} {latPutRight}], {'not','not','not','not'}, 'dud');


parInterhemiBool=bsc_extractStreamIndByName(categoryPrior,'parietal_to_parietal_interHemi');
parietalCCBool=parInterhemiBool&~forcepsMinorBool&~forcepsMajorBool&parietalCCBool&parietalEndpointsBool;

classificationOut=bsc_concatClassificationCriteria(classificationOut,'parietalCC',parietalCCBool);

frontoSeparate=bsc_planeFromROI_v2(118+[11000 12000], 'anterior',atlas);
posteriorFrontoCCLimit=bsc_planeFromROI_v2(255, 'posterior',atlas);
middleFrontalNot=bsc_planeFromROI_v2(254, 'inferior',atlas);
middleFrontalLimit=bsc_planeFromROI_v2(254, 'anterior',atlas);
anterioFrotalNot=bsc_modifyROI_v2(atlas,middleFrontalNot, middleFrontalLimit, 'posterior');

[anterioFrontalCC, anterioFrontalBool]=wma_SegmentFascicleFromConnectome(wbfg, [ {latPutLeft} {latPutRight} {frontoSeparate} {anterioFrotalNot} {posteriorFrontoCCLimit}], {'not','not','and','not','not'}, 'dud');


[middleFrontalCC, middleFrontalBool]=wma_SegmentFascicleFromConnectome(wbfg, [ {latPutLeft} {latPutRight} {frontoSeparate} {middleFrontalNot}], {'not','not','not','not'}, 'dud');

middleFrontalEndpointBool=bsc_applyEndpointCriteria(wbfg, frontoSeparate,'anterior','both');


middleFrontalInterhemiBool=bsc_extractStreamIndByName(categoryPrior,'frontal_to_frontal_interHemi');


middleFrontalBool=middleFrontalInterhemiBool&~forcepsMinorBool&~forcepsMajorBool&~parietalCCBool&middleFrontalBool&~middleFrontalEndpointBool;
anterioFrontalBool=middleFrontalInterhemiBool&~forcepsMinorBool&~forcepsMajorBool&~parietalCCBool&~middleFrontalBool&anterioFrontalBool&middleFrontalEndpointBool;
%anterioFrontalCC.fibers=wbfg.fibers(anterioFrontalBool);
%middleFrontalCC.fibers=wbfg.fibers(middleFrontalBool);

classificationOut=bsc_concatClassificationCriteria(classificationOut,'middleFrontalCC',middleFrontalBool);

if experimentalBool
classificationOut=bsc_concatClassificationCriteria(classificationOut,'anterioFrontalCC',anterioFrontalBool);
end
end