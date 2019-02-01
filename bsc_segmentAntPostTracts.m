function [classificationOut] =bsc_segmentAntPostTracts(wbfg, fsDir,varargin)
% [classificationOut] =bsc_segmentArc_Cingulum(wbfg, fsDir)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

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

categoryPrior=varargin{1};

allStreams=wbfg.fibers;
% obtain midpoints
for iFibers=1:length(allStreams)
    fiberNodeNum=round(length(allStreams{iFibers})/2);
    curStreamline=allStreams{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
end

[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);


%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

rightbool=midpoints(:,1)>0;

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

lentiLut=[12 13; 51 52];
palLut=[13;52];
thalLut=[10;49];
ventricleLut=[4;43];
wmLut=[2;41];
DCLut=[28;60];
hippLut=[17;53];

subcort=[10 12 13 17 18; 49 51 52 53 54];

    interHemiNot=bsc_makePlanarROI(atlasPath,0, 'x');
    
        [classificationOut] =bsc_segmentCingulum(wbfg, fsDir,categoryPrior);
    cingulumBool=or(classificationOut.index==find(strcmp(classificationOut.names,'rightcingulum')),classificationOut.index==find(strcmp(classificationOut.names,'leftcingulum')));
    

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% occipito roi
    
    %Planes for subcortical structures
    [thalBottom]= bsc_planeFromROI_v2(thalLut(leftright), 'inferior',atlasPath);
    palPost=bsc_planeFromROI_v2(lentiLut(leftright), 'posterior',atlasPath);
    
    %Plane for preventing crosshemispheric fiberes

    
    %planes for uncinate, also used in arcuate segmentation
    olfROI= bsc_planeFromROI_v2(164+sidenum,'inferior',atlasPath);
    antCircInsPlane=bsc_planeFromROI_v2(148+sidenum,'anterior',atlasPath);
    SupAntCirCInsPlane=bsc_modifyROI_v2(atlasPath,antCircInsPlane, olfROI, 'superior');
    PostInfOlfPlane=bsc_modifyROI_v2(atlasPath,olfROI,antCircInsPlane , 'posterior');
    
    %segment uncinate using planes
    [Unc, UncSegBool]=wma_SegmentFascicleFromConnectome(wbfg, [{SupAntCirCInsPlane},{PostInfOlfPlane},{interHemiNot}], {'and','and','not'}, 'dud');
    
    %apply midpoint and endpoint criteria
    uncMidpointBool=bsc_applyMidpointCriteria(midpoints,thalBottom,'inferior',palPost,'anterior');
    uncEndpointBool=bsc_applyEndpointCriteria(wbfg, thalBottom,'inferior','both',palPost,'anterior','both');
    
    %apply boolean criteria
    if leftright==1
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Uncinate'),UncSegBool,uncEndpointBool,uncMidpointBool,~rightbool,[efficiencyRat>.25]',categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_temporal')));
    else
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Uncinate'),UncSegBool,uncEndpointBool,uncMidpointBool,rightbool,[efficiencyRat>.25]',categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_temporal')));
    end
    %UNCINATE DONE ========================================================
    %Create anatomical Rois
    spineROI=bsc_roiFromAtlasNums(atlasPath,16,0);
    subcortROI=bsc_roiFromAtlasNums(atlasPath,subcort(leftright,:),0);
    
    %Generate anatomically defined planes
    [palTop]= bsc_planeFromROI_v2(lentiLut(leftright), 'superior',atlasPath);
    [antMargPlane] = bsc_planeFromROI_v2(147+sidenum, 'anterior',atlasPath);
    [antPalLimit] = bsc_planeFromROI_v2(palLut(leftright), 'anterior',atlasPath);
    [postPalLimit] = bsc_planeFromROI_v2(palLut(leftright), 'posterior',atlasPath);
    antParieto=bsc_planeFromROI_v2(166+sidenum,'anterior',atlasPath);
    postOrbital=bsc_planeFromROI_v2(124+sidenum,'posterior',atlasPath);
    infBrainLimit=bsc_planeFromROI_v2(121+sidenum,'inferior',atlasPath);
    
    %create modified planes as needed
    [supAntMargPlane]=bsc_modifyROI_v2(atlasPath,antMargPlane, palTop, 'superior');
    [supAntPalLimit] =bsc_modifyROI_v2(atlasPath,antPalLimit , palTop, 'superior');
    [supPostPalLimit]=bsc_modifyROI_v2(atlasPath,postPalLimit, palTop, 'superior');
    [infPostPalLimit]=bsc_modifyROI_v2(atlasPath,postPalLimit, palTop, 'inferior');
    [infAntMargPlane]=bsc_modifyROI_v2(atlasPath,antMargPlane, palTop, 'inferior');
    supOlfPlane=bsc_modifyROI_v2(atlasPath,postOrbital, olfROI, 'superior');
    infOlfPlane=bsc_modifyROI_v2(atlasPath,postOrbital, olfROI, 'inferior');
    
    %segment Using Planes
    [~, IFOFBool]=wma_SegmentFascicleFromConnectome(wbfg, [{supAntMargPlane},{supPostPalLimit},{interHemiNot},{SupAntCirCInsPlane},{antParieto},{subcortROI},{infBrainLimit}], {'not','not','not','and','and','not','not'}, 'dud');
    
    %apply boolean criteria
    if leftright==1
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'IFOF'),IFOFBool,~rightbool,[efficiencyRat>.25]',categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_occipital')));
    else
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'IFOF'),IFOFBool,rightbool,[efficiencyRat>.25]',categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_occipital')));
    end
    
    %IFOF DONE ========================================================
    %Create anatomical Rois
    insNot=bsc_roiFromAtlasNums(atlasPath,118+sidenum, 5);
    
    stemInflateA=bsc_roiFromAtlasNums(atlasPath,[168 157]+sidenum, 21);
    stemInflateB=bsc_roiFromAtlasNums(atlasPath,[141]+sidenum, 21);
    
    
    %Generate anatomically defined planes
    [thaltPost] = bsc_planeFromROI_v2(thalLut(leftright), 'posterior',atlasPath);
    [thaltTop]  = bsc_planeFromROI_v2(thalLut(leftright), 'superior',atlasPath);
    [thalBottom]= bsc_planeFromROI_v2(thalLut(leftright), 'inferior',atlasPath);
    palPost     = bsc_planeFromROI_v2(lentiLut(leftright), 'posterior',atlasPath);
    [thaltTop]  = bsc_planeFromROI_v2(thalLut(leftright), 'superior',atlasPath);
    
    %create modified planes as needed
    [infAntMargPlane] =bsc_modifyROI_v2(atlasPath,antMargPlane, palTop, 'inferior');
    [infAntPalLimit]  =bsc_modifyROI_v2(atlasPath,antPalLimit, palTop, 'inferior');
    SupAntCirCInsPlane=bsc_modifyROI_v2(atlasPath,antCircInsPlane, olfROI, 'superior');
    stemIntersect=bsc_intersectROIs(stemInflateA,stemInflateB);
    
    %segment Using Planes
    [frontalStreams, frontalStreamsBool]=wma_SegmentFascicleFromConnectome(wbfg, [{supAntMargPlane},{supAntPalLimit},{supPostPalLimit},{infAntMargPlane},{interHemiNot},{stemIntersect},{subcortROI},{spineROI},{insNot}], {'and','and','and','and','not','and','not','not','not'}, 'dud');
    
    %apply boolean criteria
    arcEndpointBool=bsc_applyEndpointCriteria(wbfg, thaltPost,'anterior','both');
    
    if leftright==1
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Arc'),frontalStreamsBool,~UncSegBool,arcEndpointBool,~rightbool,~IFOFBool,[efficiencyRat>.25]',categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_temporal')));
    else
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'Arc'),frontalStreamsBool,~UncSegBool,arcEndpointBool,rightbool,~IFOFBool,[efficiencyRat>.25]',categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_temporal')));
    end
    %Arcuate segmentation complete========================================

    [SFL1Intersection] = bsc_MultiIntersectROIs(atlasPath,13,116+sidenum, 155+sidenum, 107+sidenum );
     [SFL2Intersection] = bsc_MultiIntersectROIs(atlasPath,13,115+sidenum, 155+sidenum, 170+sidenum );
    [SFL3Intersection] = bsc_MultiIntersectROIs(atlasPath,13,112+sidenum, 150+sidenum, 153+sidenum );
    [SFL3NotIntersection] = bsc_MultiIntersectROIs(atlasPath,17,153+sidenum, 150+sidenum , 169+sidenum);
    [slf3NotTopLimit] = bsc_planeFromROI_v2(SFL3NotIntersection, 'superior',atlasPath);
    [slf3NotPostLimit] = bsc_planeFromROI_v2(SFL3NotIntersection, 'posterior',atlasPath);
    [slf3TopCut] =bsc_modifyROI_v2(atlasPath,slf3NotTopLimit, slf3NotPostLimit, 'anterior');
    
        [slf3infLimit] = bsc_planeFromROI_v2(DCLut(leftright), 'superior',atlasPath);
    [slf3infPostLimit] = bsc_planeFromROI_v2(141+sidenum, 'anterior',atlasPath);
    [slf3PostCut] =bsc_modifyROI_v2(atlasPath,slf3infLimit, slf3infPostLimit, 'posterior');
    
        
        [slf12antCutLimit] = bsc_planeFromROI_v2(155+sidenum, 'posterior',atlasPath);
    [slf12infCutLimit] = bsc_planeFromROI_v2(147+sidenum, 'inferior',atlasPath);
    [slf12exclude] =bsc_modifyROI_v2(atlasPath,slf12antCutLimit, slf12infCutLimit, 'inferior');
    
    
    cingExclude=bsc_roiFromAtlasNums(atlasPath,[108 107]+sidenum, 1);
    SL1ROIs=bsc_roiFromAtlasNums(atlasPath,[116 155]+sidenum, 1);
    SL2ROIs=bsc_roiFromAtlasNums(atlasPath,[115 155]+sidenum, 1);
    SL3ROIs=bsc_roiFromAtlasNums(atlasPath,[153 112 150]+sidenum, 1);
    periC=   bsc_roiFromAtlasNums(atlasPath,[167]+sidenum, 1);
    
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
    LentiNot=bsc_roiFromAtlasNums(atlasPath,lentiLut(leftright), 1);
    
    
    [SLF1, SLF1Bool]=wma_SegmentFascicleFromConnectome(wbfg, [{cingExclude},{SFL1Intersection},{CCintExclude},{CCvertExcludeCut}, {SFL2Intersection}, {SFL3Intersection}, {slf12exclude}], {'not','and','not','not','not','not','not'}, 'dud');
    [SLF2, SLF2Bool]=wma_SegmentFascicleFromConnectome(wbfg, [{cingExclude},{SFL2Intersection} {CCintExclude},{CCvertExcludeCut},{SFL1Intersection}, {SFL3Intersection},  {slf12exclude}], {'not','and','not','not','not','not','not'}, 'dud');
    [SLF3, SLF3Bool]=wma_SegmentFascicleFromConnectome(wbfg, [{cingExclude},{SFL3Intersection},{SFL1Intersection}, {SFL2Intersection}, {LentiNot},{slf3TopCut}], {'not','and','not','not','not','not'}, 'dud');
    SLF1Bool=SLF1Bool&categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_parietal'))&~cingulumBool;
    SLF2Bool=SLF2Bool&categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_parietal'))&~cingulumBool;
    SLF3Bool=SLF3Bool&categoryPrior.index==find(strcmp(categoryPrior.names,'frontal_to_parietal'));
    
    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF1'),SLF1Bool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF2'),SLF2Bool);
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'SLF3'),SLF3Bool);
    
    
    
end


end