function [classificationOut] =bsc_opticRadiationSeg_V7(wbfg, fsDir, expandSegBool,varargin)
% [classification] =bsc_opticRadiationSeg_V6(wbfg, fsDir, expandSegBool)
%
% This function automatedly segments components of the optic readiation
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -expandSegBool: boolean value indicating whether to perform segmentation
% of experimental/unconfirmed tracts.
%
% Outputs:
%  classificationOut:  standardly constructed classification structure
%
% (C) Daniel Bullock, 2017, Indiana University

%% parameter note & initialization
%create left/right lables
sideLabel={'left','right'};

% obtain midpoints
allStreams=wbfg.fibers;
% obtain midpoints

categoryPrior=varargin{1};

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);



atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

lentiLut=[12 13; 51 52];
palLut=[13;52];
thalLut=[10;49];
ventricleLut=[4;43];
wmLut=[2;41];
hippoLUT=[17;53];
choroLUT=[31;63];

subcort=[10 12 13 17 18; 49 51 52 53 54];
    

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    
    
    % Begin Segmentation
    %%
    % Begin segmentation of meyer's loop===================================
    % GENERATE INITIAL ANATOMICAL ROIS
    wm=bsc_roiFromAtlasNums(atlasPath,wmLut(leftright), 1);
    [thalamicROI] =bsc_roiFromAtlasNums(atlasPath,thalLut(leftright),7);
    
    %CREATE ANATOMICALLY BASED PLANES
    %postVentPlane= bsc_planeFromROI_v2([160]+sidenum,'anterior',atlasPath);
    preOccipPlane= bsc_planeFromROI_v2(172+sidenum,'posterior',atlasPath);
    %cingMargPlane= bsc_planeFromROI_v2(147+sidenum,'anterior',atlasPath);
    %thalTop      = bsc_planeFromROI_v2(thalLut(leftright),'superior',atlasPath);
    thalAnt      = bsc_planeFromROI_v2(thalLut(leftright),'anterior',atlasPath);
    thalLat      = bsc_planeFromROI_v2(thalLut(leftright),'lateral',atlasPath);
    hippTop      = bsc_planeFromROI_v2(hippoLUT(leftright),'superior',atlasPath);
     thalPost      = bsc_planeFromROI_v2(thalLut(leftright),'posterior',atlasPath);
    interHemiNot=bsc_makePlanarROI(atlasPath,0, 'x');
    
    %CREATE ROIS FROM INTERSECTIONS
    %posterior stem of Meyer
    %stemInflateA=bsc_roiFromAtlasNums(atlasPath,[174]+sidenum, 21);
    %wmIntersectA=bsc_intersectROIs(stemInflateA,wm);
    %stemInflateB=bsc_roiFromAtlasNums(atlasPath,[166]+sidenum, 21);
    %wmIntersectB=bsc_intersectROIs(stemInflateB,wm);
    %opticStem   =bsc_intersectROIs(wmIntersectB,wmIntersectA);
    
    %posterior exclusion wm
    postCingGyROI=bsc_roiFromAtlasNums(atlasPath,[110]+sidenum, 15);
    medThalWM    =bsc_intersectROIs(wm,postCingGyROI);
    postCCPlane  =bsc_planeFromROI_v2(251,'posterior',atlasPath);
    medThalWMCut =bsc_modifyROI_v2(atlasPath,medThalWM,postCCPlane , 'posterior');
      hippClip =bsc_modifyROI_v2(atlasPath,thalPost,hippTop , 'superior');
      thalLatPlost=bsc_modifyROI_v2(atlasPath,thalLat,thalPost , 'anterior');
    
    %medialPaladium exclusion wm
    palROI=bsc_roiFromAtlasNums(atlasPath,palLut(leftright), 5);
    palWM =bsc_intersectROIs(wm,palROI);
    
    %MODIFY IF NECESSARY
    %opticStemHalf=bsc_modifyROI_v2(atlasPath,opticStem,postVentPlane , 'posterior');
    %postSuperiorThalPlane=bsc_modifyROI_v2(atlasPath,cingMargPlane,thalTop , 'superior');
    
    %APPLY ENDPOINT CRITERIA
    opticBool=bsc_applyEndpointCriteria(wbfg,preOccipPlane,'posterior','one');
       %meyerEndBool=bsc_applyEndpointCriteria(wbfg,hippClip,'posterior','one');
    
       % opticMidBool=bsc_applyMidpointCriteria(wbfg,thalPost,'posterior');
    
    %SEGMENT
    [~, thalBool] = wma_SegmentFascicleFromConnectome(wbfg, [{thalamicROI}   ], {'endpoints' }, 'dud');
     [~, clipInd] = wma_SegmentFascicleFromConnectome(wbfg, [{hippClip}   ], {'not' }, 'dud');
     occipSub=bsc_extractStreamIndByName(categoryPrior,sideLabel{leftright},'occipital_to_subcortical');
     
       %[~, meyerTopInd] = wma_SegmentFascicleFromConnectome(wbfg, [{hippClip}   ], {'and' }, 'dud');
       
         [~, latThalInd] = wma_SegmentFascicleFromConnectome(wbfg, [{thalLatPlost}   ], {'and' }, 'dud');
    
    %[~, meyerBool] = wma_SegmentFascicleFromConnectome(wbfg, [{preOccipPlane} {thalamicROI} {postSuperiorThalPlane} {opticStemHalf} {thalAnt} {medThalWMCut} {palWM} {interHemiNot}], {'and','endpoints','not','and','not','not','not','not' }, 'dud');
    meyerBool=thalBool&clipInd&occipSub&latThalInd;
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'meyer'),meyerBool);

    fprintf('\n Meyers loop segmented')
   %Meyer segmentation complete============================================
   %%
   % begin Baum Segmentation
   %GENERATE ANATOMICAL ROIS
   
    %DEFINE INTERSECTION ROIS
   % baumstemInflateA=bsc_roiFromAtlasNums(atlasPath,[172]+sidenum, 21);
    %baumwmIntersectA=bsc_intersectROIs(baumstemInflateA,wm);
    %baumStemInflateB=bsc_roiFromAtlasNums(atlasPath,[166]+sidenum, 21);
    %baumWmIntersectB=bsc_intersectROIs(baumStemInflateB,wm);
    %baumStem   =bsc_intersectROIs(baumWmIntersectB,baumwmIntersectA);

    %[~, baumBool] = wma_SegmentFascicleFromConnectome(wbfg, [{preOccipPlane} {thalamicROI}  {interHemiNot} {thalAnt} {baumStem} {opticStemHalf} {medThalWMCut} {palWM} ], {'and','endpoints','not','not','and', 'not','not','not' }, 'dud');
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'baum'),thalBool,~meyerBool,occipSub,latThalInd);
    baumBool=thalBool&~meyerBool&occipSub&latThalInd;
    
    %classificationOut=bsc_reconcileClassifications(classificationOut,baumClassification);
    fprintf('\n baums loop segmented')
       %baum segmentation Complete=========================================
    %%
    %begin segmentaton of experimental tracts
    if expandSegBool
        %define relevant planes or ROIS
        %obtaining wm proximate to calcerine sulcus
        ventExpand=bsc_roiFromAtlasNums(atlasPath,ventricleLut(leftright), 5);
        ventExpandWM=bsc_intersectROIs(ventExpand,wm);
        calcerineExpand=bsc_roiFromAtlasNums(atlasPath,[145]+sidenum, 5);
        calcerineExpandWM=bsc_intersectROIs(ventExpand,calcerineExpand);
        clacerineVentWM=bsc_intersectROIs(ventExpandWM,calcerineExpandWM);
        
        hippoPlane=    bsc_planeFromROI_v2(hippoLUT(leftright),'superior',atlasPath);
        clacerineVentWMCut=bsc_modifyROI_v2(atlasPath,clacerineVentWM,hippoPlane , 'inferior');
        
 
        [~, medThalOpticBool] = wma_SegmentFascicleFromConnectome(wbfg, [{preOccipPlane} {thalamicROI} {thalAnt} {medThalWMCut} {palWM} {interHemiNot} {clacerineVentWMCut}], {'and','endpoints','not','and','not','not','not' }, 'dud');
        
        [~, calcerineOpticBool] = wma_SegmentFascicleFromConnectome(wbfg, [{preOccipPlane} {thalamicROI} {thalAnt} {medThalWMCut} {palWM} {interHemiNot} {clacerineVentWMCut} ], {'and','endpoints','not','and','not','not','and' }, 'dud');
        
        [~, palWMOpticBool] = wma_SegmentFascicleFromConnectome(wbfg, [{preOccipPlane} {thalamicROI}  {thalAnt} {medThalWMCut} {palWM} {interHemiNot}], {'and','endpoints','not','not','and','not' }, 'dud');
        
        %create classification structures
        classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'medThalOptic'),~meyerBool,opticBool,~baumBool,medThalOpticBool,~palWMOpticBool,~calcerineOpticBool);
        %classificationOut=bsc_reconcileClassifications(classificationOut,medThalOpticClassification);
        
        calcerineOpticClassification=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'calcerineOptic'),~meyerBool,opticBool,~baumBool,~medThalOpticBool,~palWMOpticBool,calcerineOpticBool);
        %classificationOut=bsc_reconcileClassifications(classificationOut,calcerineOpticClassification);
        
        palWMClassification=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'medPalOptic'),~meyerBool,opticBool,~baumBool,~medThalOpticBool,palWMOpticBool,~calcerineOpticBool);
        %classificationOut=bsc_reconcileClassifications(classificationOut,palWMClassification);
        
        fprintf('\n experimental tracts segmented')
    end
end

end