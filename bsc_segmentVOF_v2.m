function [classificationOut] =bsc_segmentVOF_v2(wbfg, fsDir,varargin)
%[classificationOut] =bsc_segmentCingulum(wbfg, fsDir,varargin)
%
% This function automatedly segments the vertical occipital fasiculus
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

categoryPrior=varargin{1};
%categoryPrior=categoryPrior{1};

%[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

wmLut=[2,41];

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');


for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    vofBottomPlane= bsc_planeFromROI_v2(162+sidenum,'superior',atlasPath);
    vofInferiorPlane= bsc_planeFromROI_v2(159+sidenum,'inferior',atlasPath);
    vofPostLimit= bsc_planeFromROI_v2(166+sidenum,'anterior',atlasPath);
    
    wmROI=bsc_roiFromAtlasNums(atlasPath,wmLut(leftright),1);
    
    vofMidpointBool=bsc_applyMidpointCriteria(wbfg,vofPostLimit,'posterior',vofBottomPlane,'superior',vofInferiorPlane,'inferior');
    [~, vofCandidateBool]=wma_SegmentFascicleFromConnectome(wbfg, [{vofBottomPlane} {vofInferiorPlane} {vofPostLimit} {wmROI}], {'and','and','not','and'}, 'dud');
        occipitalOccipitalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'occipital_to_occipital'),categoryPrior.names)))';
    
    VofBool=vofMidpointBool&vofCandidateBool&occipitalOccipitalBool;
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'VOF'),VofBool);
  
end
end