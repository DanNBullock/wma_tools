function [classificationOut] =bsc_segmentCingulum(wbfg, fsDir,varargin)
%[classificationOut] =bsc_segmentCingulum(wbfg, fsDir,varargin)
%
% This function automatedly segments the cingulum
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

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);


atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

%iterates through left and right sides
for leftright= [1,2]

    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    periCSupLim=bsc_planeFromROI_v2([156]+sidenum, 'superior',atlasPath);
    periCLatLim=bsc_planeFromROI_v2([167]+sidenum, 'lateral',atlasPath);
    [ccMidAnt] = bsc_planeFromROI_v2(252, 'anterior',atlasPath);
    
    %[cingExcPlane]=bsc_modifyROI_v2(atlasPath,ccMidAnt,periCLatLim, 'medial');
    [cingIncPlane]=bsc_modifyROI_v2(atlasPath,ccMidAnt,periCLatLim, 'lateral');
    [cingSupExcPlane]=bsc_modifyROI_v2(atlasPath,ccMidAnt,periCSupLim, 'superior');

    %goal state
    [~, cingSegBool]=wma_SegmentFascicleFromConnectome(wbfg, [{cingIncPlane}, {cingSupExcPlane} ], {'and','not'}, 'dud');
    parietoFrontalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_parietal'),categoryPrior.names)))';
    frontoTemporalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_temporal'),categoryPrior.names)))';

    classificationBool=or(parietoFrontalBool,frontoTemporalBool);
    cingRemBool=cingSegBool&classificationBool;
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'cingulum'),cingRemBool);
    
end
end

    