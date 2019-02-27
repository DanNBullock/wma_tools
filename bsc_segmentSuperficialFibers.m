function [classificationOut] =bsc_segmentSuperficialFibers(wbfg, fsDir)
%[classificationOut] =bsc_segmentCingulum(wbfg, fsDir,varargin)
%
% This function automatedly segments the supercficial fibers.
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.  Subsections may come later.

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

%[categoryPrior] =bsc_streamlineCategoryPriors_v4(wbfg, fsDir,2)

%categoryPrior=categoryPrior{1};

%[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

wmLut=[2,41];

streamLengthLimit=30;

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

[inflatedAtlas] =bsc_inflateLabels(fsDir,2);

greyMatterROIS=[[101:1:175]+12000 [101:1:175]+11000];

[greyROI] =bsc_roiFromAtlasNums(inflatedAtlas,greyMatterROIS,1);

midpointBool=bsc_midpointROISegment(wbfg,greyROI);


for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    wmROIDegrade=bsc_roiFromAtlasNums(inflatedAtlas,wmLut(leftright),1);
    wmROIFull=bsc_roiFromAtlasNums(atlasPath,wmLut(leftright),1);
    
    if leftright==1
        OtherwmROIFull=bsc_roiFromAtlasNums(atlasPath,wmLut(2),5);
    else
        OtherwmROIFull=bsc_roiFromAtlasNums(atlasPath,wmLut(1),5);
    end
    z = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), wbfg.fibers, 'UniformOutput', true);
    
    [superficialFibersDegrade, superficialFibersDegradeBool]=wma_SegmentFascicleFromConnectome(wbfg, [{wmROIDegrade} {OtherwmROIFull}], {'not','not'}, 'dud');
    %[superficialFibersFull, superficialFibersFullBool]=wma_SegmentFascicleFromConnectome(wbfg, [{wmROIFull} {OtherwmROIFull}], {'not','not'}, 'dud');
    
    
    %selectFibers.fibers=wbfg.fibers(midpointBool&superficialFibersDegradeBool);
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Ufibers'),z<streamLengthLimit,midpointBool,superficialFibersDegradeBool);
    
end
end