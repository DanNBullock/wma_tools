function [classificationOut] =bsc_segmentCingulum_v3(wbfg,atlasPath,varargin)
%
% This function automatedly segments the cingulum
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
thalIDs=[10,49];

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

%iterates through left and right sides
for leftright= [1,2]
    sidenum=10000+leftright*1000;
    thalTop=bsc_planeFromROI_v2(thalIDs(leftright), 'superior',atlasPath);
    ccPostLimit=bsc_planeFromROI_v2(251, 'posterior',atlasPath);
    ccAntLimit=bsc_planeFromROI_v2(254, 'anterior',atlasPath);
    
    %carve out the area around and above the cc
    [postCCtopThal]=bsc_modifyROI_v2(atlasPath,ccPostLimit, thalTop, 'superior');
    [ccInterior1]=bsc_modifyROI_v2(atlasPath,thalTop, ccPostLimit, 'anterior');
    [ccInterior2]=bsc_modifyROI_v2(atlasPath,ccInterior1, ccAntLimit, 'posterior');

    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    
    periCSupLim=bsc_planeFromROI_v2([156]+sidenum, 'superior',atlasPath);
    periCLatLim=bsc_planeFromROI_v2([167]+sidenum, 'lateral',atlasPath);
    [ccMidAnt] = bsc_planeFromROI_v2(252, 'anterior',atlasPath);
    
    
    [cingIncPlane]=bsc_modifyROI_v2(atlasPath,ccMidAnt,periCLatLim, 'medial');
    
    %[cingExcPlane]=bsc_modifyROI_v2(atlasPath,ccMidAnt,periCLatLim, 'medial');
    
    [cingSupExcPlane]=bsc_modifyROI_v2(atlasPath,ccMidAnt,periCSupLim, 'superior');
    
    %goal state
    [~, cingSegBool]=wma_SegmentFascicleFromConnectome(wbfg, [{cingIncPlane}, {cingSupExcPlane}, {ccInterior2} ], {'and','not','not'}, 'dud');
    
    parietoFrontalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel(leftright),'frontal_to_parietal'));
    frontoTemporalBool=bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel(leftright),'frontal_to_temporal'));
    
    classificationBool=or(parietoFrontalBool,frontoTemporalBool);
    cingRemBool=cingSegBool&classificationBool;
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'cingulum'),cingRemBool);
    
end
end

    
