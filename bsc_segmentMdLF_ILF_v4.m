function [classificationOut] =bsc_segmentMdLF_ILF_v4(wbfg,atlas,categoryPrior)
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
% (C) Daniel Bullock, 201*, Indiana University

%% parameter note & initialization


% caudate ROI specification
caudateLUT=[11,50];
amygdlaIDs=[18,54];
thalIDs=[10,49];
subCortIDs=[5 12 13 17 18 ;44 51 52 53 54];

%create left/right lables
sideLabel={'left','right'};

%atlas=fullfile(fsDir,'mri/aparc.a2009s+aseg.nii.gz');

% obtain midpoints
allStreams=wbfg.fibers;

%initialize classification structure
classification=[];
classification.names=[];
classification.index=zeros(length(wbfg.fibers),1);
classificationILF=classification;
classificationMDLFang=classification;
classificationMDLFspl=classification;

[inflatedAtlas] =bsc_inflateLabels(fsDir,2);

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% occipito roi
    %generates the roi for the occipito-parietal regions corresponding to
    %MdLF --no 145
    [mergedOCTROI] =bsc_roiFromAtlasNums(inflatedAtlas,[120,111,166, 143,119,158,122,102,119,159]+sidenum,1);
    
    %curtail somewhat; provide anterior limit, here we use anterior of IPS
    %as limit
    [antLimit] = bsc_planeFromROI_v2(157+sidenum, 'anterior',atlas);
    [mergedOCTROI]=bsc_modifyROI_v2(atlas,mergedOCTROI, antLimit, 'posterior');
    
    %% parietal roi
    %creates ROI just for the Parietal area
    [mergedParietalROI] =bsc_roiFromAtlasNums(inflatedAtlas,[157, 127, 168, 136, 126, 125]+sidenum,5);
    
    %% lateral temporal roi
    %generates the roi for the lateral-temporal regions corresponding to
    %MdLF
    [mergedLatTempROI] =bsc_roiFromAtlasNums(inflatedAtlas,[134, 144, 174,135]+sidenum,9);
    
    %impliments a cutoff to ensure that the temporal roi coordinates are
    %anterior of the y=-15 plane
    %[amygdlaROI] =bsc_roiFromAtlasNums(atlas,[134, 122, 144, 174, 114, 135]+sidenum,1,9);
    
    %use amygdala to institute planar limit on mergedLatTempROI
    [amygdalaPost] = bsc_planeFromROI_v2(amygdlaIDs(leftright), 'posterior',atlas);
    [mergedLatTempROI]=bsc_modifyROI_v2(atlas,mergedLatTempROI, amygdalaPost, 'anterior');
    
    %% begin segmenting
    %bsc_extractStreamIndByName(classification,tractName)
    occipitoTemporalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'occipital_to_temporal'));
    
    %get rid of ILF streamlines above anterior temporal area
    [lingGyAnt] = bsc_planeFromROI_v2(162+sidenum, 'anterior',atlas);
    [posteriorILFAnt] = bsc_planeFromROI_v2(110+sidenum, 'posterior',atlas);
    
    ILFAntEndpointBool=bsc_applyEndpointCriteria(wbfg,lingGyAnt,'anterior','one');
    ILFPostEndpointBool=bsc_applyEndpointCriteria(wbfg,posteriorILFAnt,'posterior','one');
    
    insPost=bsc_planeFromROI_v2(150+sidenum, 'posterior',atlas);
    tempTransVTop=bsc_planeFromROI_v2(133+sidenum, 'superior',atlas);
    TopArcAnd=bsc_modifyROI_v2(atlas,insPost, tempTransVTop, 'superior');
    
    cingRemoveROI=bsc_roiFromAtlasNums(inflatedAtlas,[109 110]+sidenum,5);
    
    
    [~, ILFBool]=wma_SegmentFascicleFromConnectome(wbfg, [{TopArcAnd} {cingRemoveROI} ], {'not', 'not'}, 'dud');
    
    
    [classificationILF]=bsc_concatClassificationCriteria(classificationILF,strcat(sideLabel{leftright},'ILF'),occipitoTemporalBool,ILFAntEndpointBool,ILFPostEndpointBool,ILFBool);
    
    
    [~, MDLFind]=bsc_tractByEndpointROIs(wbfg, [{mergedParietalROI},{mergedLatTempROI}]);
    
    %Make ROI for subcortical areas, later to be used to exclude
    %streamlines passing through these areas
    [subCortROI] = bsc_roiFromAtlasNums(atlas,subCortIDs(leftright,:),1);
    
    [anteriorMDLFLimit] = bsc_planeFromROI_v2(172+sidenum, 'anterior',atlas);
    [latFis] = bsc_planeFromROI_v2(150+sidenum, 'superior',atlas);
    [mdlfFrontRm]=bsc_modifyROI_v2(atlas,anteriorMDLFLimit, latFis, 'superior');
    
    [~, removeMDLFIND]=wma_SegmentFascicleFromConnectome(wbfg, [{subCortROI},{mdlfFrontRm}], {'not','not'}, 'dud');

    temporoParietalBool= bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'parietal_to_temporal'))';
    
    %Find WM for spl
    %extreme inflation of precuneus
    splMid=bsc_roiFromAtlasNums(atlas,130+sidenum,27);
    %extreme inflation of IPS
    splLat=bsc_roiFromAtlasNums(atlas,157+sidenum,27);
    %find out where they intersect to get white matter of superior parietal
    %lobule
    splWM=bsc_intersectROIs(splMid,splLat);
    
    %find streamlines that intersect the superior parietal lobule white
    %matter
    [~, splIND]=wma_SegmentFascicleFromConnectome(wbfg, [{splWM}], {'and'}, 'dud');
    
    %midpoints should be below the top of the thalamus if they are posterior to the
    %amygdala and simply should not be anterior of the amygdala in any
    %event.
    subcortTop=bsc_planeFromROI_v2(subCortROI,'superior',atlas) ;
    [badMidpointsA]=bsc_applyMidpointCriteria(wbfg, anteriorMDLFLimit,'anterior',subcortTop,'superior',amygdalaPost,'posterior');
    [badMidpointsB]=bsc_applyMidpointCriteria(wbfg,amygdalaPost,'anterior');
    
    %apply the criteria to find the streamlines corresponding to the
    %MDLFang and MDLFspl.  Note that the main thing that is distinguishing
    %them in this segmentation is the
    [classificationMDLFang]=bsc_concatClassificationCriteria(classificationMDLFang,strcat(sideLabel{leftright},'MDLFang'),MDLFind,removeMDLFIND,~badMidpointsA,~badMidpointsB,~splIND,temporoParietalBool);
    [classificationMDLFspl]=bsc_concatClassificationCriteria(classificationMDLFspl,strcat(sideLabel{leftright},'MDLFspl'),MDLFind,removeMDLFIND,~badMidpointsA,~badMidpointsB,splIND,temporoParietalBool);
    
end

%Compile classifications into one single structure
classificationOut=bsc_reconcileClassifications(classification,classificationILF);
classificationOut=bsc_reconcileClassifications(classificationOut,classificationMDLFang);
classificationOut=bsc_reconcileClassifications(classificationOut,classificationMDLFspl);
end