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

%[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);


atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

lentiLut=[12 13; 51 52];
palLut=[13;52];
wmLut=[2;41];

%iterates through left and right sides
for leftright= [1,2]
    
    [CCintExclude] = bsc_planeFromROI_v2(254, 'inferior',atlasPath);
    [CCvertExclude] = bsc_planeFromROI_v2(252, 'posterior',atlasPath);
    [CCvertExcludeLimit] = bsc_planeFromROI_v2(254, 'superior',atlasPath);
    [CCvertAnteriorExcludeLimit] = bsc_planeFromROI_v2(253, 'posterior',atlasPath);
    [CCvertExcludeCut] =bsc_modifyROI_v2(atlasPath,CCvertExclude, CCvertExcludeLimit, 'inferior');
    [CCvertAnteriorExcludeCut] =bsc_modifyROI_v2(atlasPath,CCvertAnteriorExcludeLimit, CCvertExcludeLimit, 'inferior');
    [CCExcludeAntLimit] = bsc_planeFromROI_v2(255, 'anterior',atlasPath);
    [CCExcludePosteriorLimit] = bsc_planeFromROI_v2(251, 'posterior',atlasPath);
    [CCintExclude] =bsc_modifyROI_v2(atlasPath,CCintExclude, CCExcludeAntLimit, 'posterior');
    [CCintExclude] =bsc_modifyROI_v2(atlasPath,CCintExclude, CCExcludePosteriorLimit, 'anterior');
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    cingLatBorder=bsc_planeFromROI_v2([157]+sidenum, 'medial',atlasPath);
    cingEndpointBool=bsc_applyEndpointCriteria(wbfg,cingLatBorder,'lateral','neither');
    
    [palTop]= bsc_planeFromROI_v2(lentiLut(leftright), 'superior',atlasPath);
    [antMargPlane] = bsc_planeFromROI_v2(147+sidenum, 'anterior',atlasPath);
    [antPalLimit] = bsc_planeFromROI_v2(palLut(leftright), 'anterior',atlasPath);
    
    %create modified planes as needed
    [supAntMargPlane]=bsc_modifyROI_v2(atlasPath,antMargPlane, palTop, 'superior');
    [supAntPalLimit] =bsc_modifyROI_v2(atlasPath,antPalLimit , palTop, 'superior');
    
    cingShortExcludeBool1=bsc_applyEndpointCriteria(wbfg,supAntMargPlane,'anterior','both');
    cingShortExcludeBool2=bsc_applyEndpointCriteria(wbfg,supAntPalLimit,'posterior','both');
    cingShortExcludeBoolBoth=or(cingShortExcludeBool1,cingShortExcludeBool2);
    
    wm=bsc_roiFromAtlasNums(atlasPath,wmLut(leftright), 1);
    

    
    cingMargWM=bsc_roiFromAtlasNums(atlasPath,[147]+sidenum, 21);
    CingMargWMIntersectA=bsc_intersectROIs(cingMargWM,wm);
    
    %goal state
    [frontalRem2, frontalStreamsRemBool2]=wma_SegmentFascicleFromConnectome(wbfg, [{CingMargWMIntersectA} {CCvertExcludeCut} {CCintExclude}], {'and','not','not'}, 'dud');
      parietoFrontalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_parietal'),categoryPrior.names)))';
       frontoTemporalBool=(categoryPrior.index==find(strcmp(strcat(sideLabel(leftright),'frontal_to_temporal'),categoryPrior.names)))';
     
      
    
    classificationBool=or(parietoFrontalBool,frontoTemporalBool);
    frontalStreamsRemBool2=~cingShortExcludeBoolBoth&cingEndpointBool&frontalStreamsRemBool2&classificationBool;
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'cingulum'),frontalStreamsRemBool2);
 
    %frontalRem2.fibers=wbfg.fibers(frontalStreamsRemBool2);
    %     figure

end
end
    
    