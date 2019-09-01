function [classificationOut] =bsc_segmentCingulum_v3(wbfg,atlas,varargin)
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

thalIDs=[10,49];

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

%atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

%iterates through left and right sides
for leftright= [1,2]
    sidenum=10000+leftright*1000;
    thalTop=bsc_planeFromROI_v2(thalIDs(leftright), 'superior',atlas);
    ccPostLimit=bsc_planeFromROI_v2(251, 'posterior',atlas);
    ccAntLimit=bsc_planeFromROI_v2(254, 'anterior',atlas);
    
    %carve out the area around and above the cc
    [postCCtopThal]=bsc_modifyROI_v2(atlas,ccPostLimit, thalTop, 'superior');
    [ccInterior1]=bsc_modifyROI_v2(atlas,thalTop, ccPostLimit, 'anterior');
    [ccInterior2]=bsc_modifyROI_v2(atlas,ccInterior1, ccAntLimit, 'posterior');

    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    
    periCSupLim=bsc_planeFromROI_v2([156]+sidenum, 'superior',atlas);
    periCLatLim=bsc_planeFromROI_v2([167]+sidenum, 'lateral',atlas);
    [ccMidAnt] = bsc_planeFromROI_v2(252, 'anterior',atlas);
    
    
    [cingIncPlane]=bsc_modifyROI_v2(atlas,ccMidAnt,periCLatLim, 'medial');
    [cingSupExcPlane]=bsc_modifyROI_v2(atlas,ccMidAnt,periCSupLim, 'superior');
    
    %goal state
    [~, cingSegBool]=wma_SegmentFascicleFromConnectome(wbfg, [{cingIncPlane}, {cingSupExcPlane}, {ccInterior2} ], {'and','not','not'}, 'dud');
    
    parietoFrontalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_parietal'));
    frontoTemporalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_temporal'));
    
    classificationBool=or(parietoFrontalBool,frontoTemporalBool);
    cingRemBool=cingSegBool&classificationBool;
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'cingulum'),cingRemBool);
    
end
end

    