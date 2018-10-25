function  [classificationOut] = bsc_segpArcTPC(wbfg, fsDir)
% [classificationOut] = bsc_segpArcTPC(wbfg, fsDir)
%
% This function automatedly segments the posterior arcuate and temporo-parietal
% connection from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.
%
% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
%
% Outputs:
% classification:  a classification structure with .name and .indexes
% fields
%
% (C) Daniel Bullock, 2018, Indiana University
%% parameter note & initialization

%these 3 digit numbers correspond to the last 3 digits of the DK 2009
%freesurfer look up table numbers.
%(see:  https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT)
%when added with sidenum you get the actual roi ID nums

parietalROIs3=[157, 127, 168, 136, 126, 125];
temporalROIs3=[121, 161, 137, 162, 138, 173];

%maybe play with this if something isn't to your liking, it corresponds to
%the smoothing kernel used for the parietal and temporal ROIs

smoothParameter=1;

%initialize classification structure
classification=[];
classification.names=[];
classification.index=zeros(length(wbfg.fibers),1);
classificationpArc=classification;
classificationTPC=classification;

%create left/right lables
sideLabel={'left','right'};

%indexes for the putamen
putInd=[12,51];
ventInd=[4,43];

%% actual segmentation

%iterates through left and right sides
for leftright= [1,2]
    %create temporary classification structure
    classificationTemp=[];
    classificationTemp.names=[];
    classificationTemp.index=zeros(length(wbfg.fibers),1);
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% parietal ROI
    %generates the roi for the parietal regions corresponding to the pArc
    %and TPC
    
    [mergedParietalROI] =bsc_roiFromFSnums(fsDir,parietalROIs3+sidenum,1,smoothParameter);
    mergedParietalROI.name='parietalROI';
    
    %% temporal ROI
    %generates the roi for the temporal regions corresponding to the pArc
    %and TPC
    [mergedTemporalROI] =bsc_roiFromFSnums(fsDir,temporalROIs3+sidenum,1,smoothParameter);
    mergedTemporalROI.name='temporalROI';
    
    %create a demarcation based on the posterior putamen
    [postPut] = bsc_planeFromROI(putInd(leftright), 'posterior',fsDir);
    
    %limit temporal ROI to those areas posterior to the putamen
    [mergedTemporalROI]=bsc_modifyROI(fsDir,mergedTemporalROI, postPut, 'posterior');
    
    %% Additional Requirements
    %Require that it pass through a plane at the bottom of the ips (156
    %isn't the ips itself but is reliably slightly lower than it
    [supPlane] = bsc_planeFromROI(156+sidenum, 'inferior',fsDir);
    
    %Require that it pass through a plane at the slightly inferior to
    %the temporoparietal junction
    [infPlane] = bsc_planeFromROI(136+sidenum, 'inferior',fsDir);
    
    %create a plane at the posterior of the frontal lobe, to prevent
    %streamlines from proceeding too far anterior
    frontPost=   bsc_planeFromROI(116+sidenum, 'posterior',fsDir);
    %ensure it only impacts the superior terminations
    frontPost=bsc_modifyROI(fsDir,frontPost, supPlane, 'superior');
    %apply it to get boolean output
    [~, vertInd] =wma_SegmentFascicleFromConnectome(wbfg, [{supPlane},{infPlane},{frontPost}], {'and','and','not'}, 'verticalTraverse');
    
    %% segment
    %find the streamlines that run between the primary ROIs
    [~, endPointInd]=  bsc_tractByEndpointROIs(wbfg, [{mergedParietalROI},{mergedTemporalROI}]);
    
    %extreme inflation of precuneus
    tpcMid=bsc_roiFromFSnums(fsDir,130+sidenum,1,27);
    %extreme inflation of IPS
    tpcLat=bsc_roiFromFSnums(fsDir,157+sidenum,1,27);
    %find out where the intersect to get white matter of superior parietal
    %lobule
    tpcWM=bsc_intersectROIs(tpcMid,tpcLat);
    
    %find top of lingual sulcus
    fusTop=bsc_planeFromROI(162+sidenum,'superior',fsDir);
    %find posterior of ventricles
    ventPost=bsc_planeFromROI(ventInd(leftright),'posterior',fsDir);
    %modify the plane at the top of the fusiform gyrus to ensure it only
    %proceeds posteriorly from the ventricle.  This is to prevent some
    %errant/spurrous tracts in the pArc
    fusTopPost=bsc_modifyROI(fsDir,fusTop, ventPost, 'posterior');
    
    %find the streamlines that intersect with the SPL
    [~, splInd] =wma_SegmentFascicleFromConnectome(wbfg, [{tpcWM}], {'and',}, 'spl');
    
    %find streamlines that intersect the fusTopPost exclusion
    [~, postRemove] =wma_SegmentFascicleFromConnectome(wbfg, [{fusTopPost}], {'not',}, 'postNot');
    
    %apply various criteria to get classification structure
    [classificationpArc]=bsc_concatClassificationCriteria(classificationpArc,strcat(sideLabel(leftright),'pArc'),endPointInd,vertInd,~splInd,postRemove);
    [classificationTPC]=bsc_concatClassificationCriteria(classificationTPC,strcat(sideLabel(leftright),'TPC'),endPointInd,vertInd,splInd);
end

%merge classification structures
[classificationOut] = bsc_reconcileClassifications(classificationpArc,classificationTPC);
end
