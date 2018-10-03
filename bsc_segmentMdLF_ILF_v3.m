function [classificationOut] =bsc_segmentMdLF_ILF_v3(wbfg, fsDir)
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

% obtain midpoints
for iFibers=1:length(wbfg.fibers)
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
end

%initialize classification structure
classification=[];
classification.names=[];
classification.index=zeros(length(wbfg.fibers),1);
classificationILF=classification;
classificationMDLFang=classification;
classificationMDLFspl=classification;

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% occipito roi
    %generates the roi for the occipito-parietal regions corresponding to
    %MdLF --no 145
    [mergedOCTROI] =bsc_roiFromFSnums(fsDir,[120,111,166, 143,119,158,122,102,119,159]+sidenum,1,5);
    
    %curtail somewhat; provide anterior limit, here we use anterior of IPS
    %as limit
    [antLimit] = bsc_planeFromROI(157+sidenum, 'anterior',fsDir);
    [mergedOCTROI]=bsc_modifyROI(fsDir,mergedOCTROI, antLimit, 'posterior');
    
    %% parietal roi
    %creates ROI just for the Parietal area
    [mergedParietalROI] =bsc_roiFromFSnums(fsDir,[157, 127, 168, 136, 126, 125]+sidenum,1,5);
    
    %% lateral temporal roi
    %generates the roi for the lateral-temporal regions corresponding to
    %MdLF
    [mergedLatTempROI] =bsc_roiFromFSnums(fsDir,[134, 144, 174,135]+sidenum,1,9);
    
    %impliments a cutoff to ensure that the temporal roi coordinates are
    %anterior of the y=-15 plane
    %[amygdlaROI] =bsc_roiFromFSnums(fsDir,[134, 122, 144, 174, 114, 135]+sidenum,1,9);
    
    %use amygdala to institute planar limit on mergedLatTempROI
    [amygdalaPost] = bsc_planeFromROI(amygdlaIDs(leftright), 'posterior',fsDir);
    [mergedLatTempROI]=bsc_modifyROI(fsDir,mergedLatTempROI, amygdalaPost, 'anterior');
    
    %% begin segmenting
    
    %segment ILF by using parieto-ocipital ROIs and lateral temporal ROIs
    [~, ILFFiberBoolVec]=bsc_tractByEndpointROIs(wbfg, [{mergedOCTROI},{mergedLatTempROI}]);
    
    %get rid of ILF streamlines above anterior temporal area 
    [tempOutRemAnt] = bsc_planeFromROI(175+sidenum, 'anterior',fsDir);
    [tempOutRemSup] = bsc_planeFromROI(175+sidenum, 'inferior',fsDir);
    [tempOutRemAntMod]=bsc_modifyROI(fsDir,tempOutRemAnt, tempOutRemSup, 'superior');
    
    %Get rid of streamlines that curve back anteriorly in the posterior ILF
    [medialExcludeA] = bsc_planeFromROI(167+sidenum, 'posterior',fsDir);
    [medialExcludeB] = bsc_planeFromROI(thalIDs(leftright), 'lateral',fsDir);
    [medialExcludeC]=bsc_modifyROI(fsDir,medialExcludeA, medialExcludeB, 'lateral');
    
    %Find indexes of streamlines that intersect the previous to planar ROIs
    [~, removeILFIND]=wma_SegmentFascicleFromConnectome(wbfg, [{tempOutRemAntMod},{medialExcludeC}], {'not','not'}, 'dud');
    
    %Apply inclusion and inclusion criteria to find indexes of ILF
    [classificationILF]=bsc_concatClassificationCriteria(classificationILF,strcat(sideLabel{leftright},'ILF'),ILFFiberBoolVec,removeILFIND);
    
    %Now we move to the MDLF
    %Do initial segmentation to get combined MDLF indexes
    [~, MDLFind]=bsc_tractByEndpointROIs(wbfg, [{mergedParietalROI},{mergedLatTempROI}]);
    
    %Make ROI for subcortical areas, later to be used to exclude
    %streamlines passing through these areas
    [subCortROI] = bsc_roiFromFSnums(fsDir,subCortIDs(leftright,:));
    
    %Create an exclusion plane somewhat similar to tempOutRemAntMod, but
    %slightly posterior for the purposes of removing streamlines that curve
    %around into frontal lobe.  These streamlines are likely part of the
    %arcuate.
    [anteriorMDLFLimit] = bsc_planeFromROI(172+sidenum, 'anterior',fsDir);
    [latFis] = bsc_planeFromROI(150+sidenum, 'superior',fsDir);
    [mdlfFrontRm]=bsc_modifyROI(fsDir,anteriorMDLFLimit, latFis, 'superior');
    
    %find streamlines that should be excluded from MDLF
    [~, removeMDLFIND]=wma_SegmentFascicleFromConnectome(wbfg, [{subCortROI},{mdlfFrontRm}], {'not','not'}, 'dud');
   
    %Find WM for spl
    %extreme inflation of precuneus
    splMid=bsc_roiFromFSnums(fsDir,130+sidenum,1,27);
    %extreme inflation of IPS
    splLat=bsc_roiFromFSnums(fsDir,157+sidenum,1,27);
    %find out where they intersect to get white matter of superior parietal
    %lobule
    splWM=bsc_intersectROIs(splMid,splLat);
    
    %find streamlines that intersect the superior parietal lobule white
    %matter
    [~, splIND]=wma_SegmentFascicleFromConnectome(wbfg, [{splWM}], {'and'}, 'dud');
    
    %midpoints should be below the top of the thalamus if they are posterior to the
    %amygdala and simply should not be anterior of the amygdala in any
    %event.
     subcortTop=bsc_planeFromROI(subCortROI,'superior',fsDir) ;
    [badMidpointsA]=bsc_applyMidpointCriteria(midpoints, anteriorMDLFLimit,'anterior',subcortTop,'superior',amygdalaPost,'posterior');
    [badMidpointsB]=bsc_applyMidpointCriteria(midpoints,amygdalaPost,'anterior');
    
    %apply the criteria to find the streamlines corresponding to the
    %MDLFang and MDLFspl.  Note that the main thing that is distinguishing
    %them in this segmentation is the
    [classificationMDLFang]=bsc_concatClassificationCriteria(classificationMDLFang,strcat(sideLabel{leftright},'MDLFang'),MDLFind,removeMDLFIND,~badMidpointsA,~badMidpointsB,~splIND);
    [classificationMDLFspl]=bsc_concatClassificationCriteria(classificationMDLFspl,strcat(sideLabel{leftright},'MDLFspl'),MDLFind,removeMDLFIND,~badMidpointsA,~badMidpointsB,splIND);
    
end

%Compile classifications into one single structure
classificationOut=bsc_reconcileClassifications(classification,classificationILF);
classificationOut=bsc_reconcileClassifications(classificationOut,classificationMDLFang);
classificationOut=bsc_reconcileClassifications(classificationOut,classificationMDLFspl);                
    end