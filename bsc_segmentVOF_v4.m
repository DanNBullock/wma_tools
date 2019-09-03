function [classificationOut] =bsc_segmentVOF_v4(wbfg,atlas,varargin)
% This function automatedly segments the vertical occipital fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.  Subsections may come later.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -varargin: priors from previous steps

% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University

sideLabel={'left','right'};
categoryPrior=varargin{1};

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

wmLut=[2,41];

inferiorOccipROInums=[143 102 122 152 161 162];
superiorOccipROInums=[120 121 158 111 166 159];
neitherROInums=[145 160];

for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    inferiorOccipROI=bsc_roiFromAtlasNums(atlas,inferiorOccipROInums+sidenum,1);
    [~, inferiorBool]=  bsc_tractByEndpointROIs(wbfg, [{inferiorOccipROI} {inferiorOccipROI}]);
    
        superiorOccipROI=bsc_roiFromAtlasNums(atlas,superiorOccipROInums+sidenum,1);
    [~, superiorBool]=  bsc_tractByEndpointROIs(wbfg, [{superiorOccipROI} {superiorOccipROI}]);
    
    [~,criteriaBool]= bsc_tractByEndpointROIs(wbfg, [{inferiorOccipROI} {superiorOccipROI}]);
    
    [botthEndpointsInf]=bsc_applyEndpointCriteria(wbfg, [0 0 0], 'inferior','both');
    
    [supEndpointsCriter]=bsc_applyEndpointCriteria(wbfg, [0 0 3], 'superior','one');
    [infEndpointsCriter]=bsc_applyEndpointCriteria(wbfg, [0 0 -3], 'inferior','one');
    
    vofBottomPlane= bsc_planeFromROI_v2(162+sidenum,'superior',atlas);
    vofInferiorPlane= bsc_planeFromROI_v2(159+sidenum,'inferior',atlas);
    vofPostLimit= bsc_planeFromROI_v2(111+sidenum,'anterior',atlas);
    
    wmROI=bsc_roiFromAtlasNums(atlas,wmLut(leftright),1);
    
    vofMidpointBool=bsc_applyMidpointCriteria(wbfg,vofPostLimit,'posterior',vofBottomPlane,'superior',vofInferiorPlane,'inferior');
    [~, vofCandidateBool]=wma_SegmentFascicleFromConnectome(wbfg, [{vofPostLimit}], {'not'}, 'dud');
    %should help fix failures to find vof
    [occipitalOccipitalBool] = or(bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'occipital_to_occipital')),bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'occipital_to_occipital_ufiber')));
    %occipitalOccipitalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'occipital_to_occipital'),categoryPrior.names)))';
    
    VofBool=~inferiorBool'&~superiorBool'&occipitalOccipitalBool&criteriaBool'&vofCandidateBool&supEndpointsCriter&infEndpointsCriter;
    VofBoolDud=~inferiorBool'&~superiorBool'&occipitalOccipitalBool&criteriaBool'&vofCandidateBool&~botthEndpointsInf;

    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'VOF'),VofBool);
  
end
end
