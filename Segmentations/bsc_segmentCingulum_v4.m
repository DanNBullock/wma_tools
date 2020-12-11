function [classificationOut] =bsc_segmentCingulum_v4(wbfg, atlas, categoryClassification)
% [classificationOut] = bsc_segmentCingulum_v4(wbfg, fsDir, categoryClassification)
%
% This is a standalone script for the segmentation of the cingulum.  Due to
% the complexity of segmenting this structure, it is given its own script
% here.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -atlas: path to the atlas  used in the segmentation 
% -categoryClassification: the classification structure resulting from the
% classification segmentation.  Done outside of this function to avoid
% doing it repeatedly

% Outputs:
% -classificationOut:  standardly constructed classification structure

% (C) Daniel Bullock, 2020, Indiana University

%% parameter notes & initialization

%create left/right labels.  For use with naming conventions later.
sideLabel={'left','right'};

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
%make sure it is the same length as the input
classificationOut.index=zeros(length(wbfg.fibers),1);

%if the passed variable for atlas is the path to the atlas, load it, if
%not, do nothing
if ischar(atlas)
atlas=niftiRead(atlas);
else
    %do nothing
end

%Set some initial rois that don't follow a good convention.  The reason we
%are doing this is that typically, when use aparc.a2009s we can designate
%left or right by adding 12000 or 11000 to a three digit number
%corresponding to a cortical roi.  Subcortical rois (which are essential to
%anatomically based segmentations) do not follow this convention, and so
%they must be done in the following way.  We can thus select between right
%and left rois by indexing into the variable with {1} or {2}, as set by the
%subsequent leftright variable
lentiLut=[12 13; 51 52];
palLut=[13;52];
thalLut=[10;49];
ventricleLut=[4;43];
wmLut=[2;41];
DCLut=[28;60];
hippLut=[17;53];
amigLut=[18;54];
nACLut=[26;58];


%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% Cingulum
    % DOI
    % The cingulum is a deep (within the brain) white matter structure that
    % runs above the corpus callosum.  Various reports have indicated that
    % there may be many subsections to this structure.  For now, this
    % semgenation will extract 2 subsections of it.  One which connects
    % parietal and posterior frontal regions to anterior frontal regions,
    % and another which curves around the anterior corpus callosum and
    % terminates in the anterior cubcallosal region.
%========================================================================= 

%1.  Cingulum - ESTABLISH CATEGORY CRITERIA    
%  Given the complex substructuring of the cingulum, it is not as though
%  there is a single category that it belongs to.  As such we have to make
%  amalgums of several categories.  This has the unfortunate side effect of
%  muddling the clarity that is typically provided by the category
%  segmentation.  As a consequences of this, the subsequent segmentation
%  steps will be more complicated.

% Regions for the fronto-parietal cingulum:
    frontoParietalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_parietal'));
    frontoFrontalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_frontal'));

% Regions for the subcallosal cingulum:
    frontoAntSubCortBool= bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_caudateNAc_to_frontal'));
    parietoAntSubCortBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_caudateNAc_to_parietal'));
%  For both:
    pericollosallBool=  bsc_extractStreamIndByName(categoryClassification,'pericallosal');
    
    % Regions for the posterior cingulum:
    occipitoHippBool= bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_hippAmig_to_occipital'));
    parietoAntSubCortBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_caudateNAc_to_parietal'));

%  applying them as a conjunct
    frontalConjunctBool=any([frontoParietalBool,frontoFrontalBool,pericollosallBool],2);
    subcallosalConjunctBool=any([frontoAntSubCortBool, parietoAntSubCortBool , pericollosallBool],2);
    
%--------------------------------------------------------------------------
%========================================================================= 

%2.  Cingulum - ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA  
% In the case of the fronto-parietal cingulum, additional endpoint criteria
% are needed in order to prevent shorter / non - cingulum streamlines from
% being segmented.  As with previous segmenations, this involves
% establishing an anterior and posterior border for these endpoints.
%--------------------------------------------------------------------------

% bordered endpoints
    
% here we use the posterior border of the relatively central 3rd ventricle
% to establish a posterior border of the endpoint criteria (we require at
% one endpoint to be posterior to this).  This anatomical landmark tends to
% be slightly anterior of the posterior border of the midbrain
    thirdVentPosterior=bsc_planeFromROI_v2(14,'posterior',atlas);
% here we do the same for the anterior endpoints, using the posterior
% boder of the subcallosal gyrus
    anteriorCurveLimit =bsc_planeFromROI_v2(132+sidenum, 'posterior',atlas);
    
    posteriorEndpointBool=bsc_applyEndpointCriteria(wbfg,thirdVentPosterior,'posterior','one');
    anteriorEndpointBool=bsc_applyEndpointCriteria(wbfg,anteriorCurveLimit,'anterior','one');
    
%--------------------------------------------------------------------------
% exclusion of superior frontal white matter
% In order to distinguish this tract from the SLF we need to exclude
% streamlines which occupy the superior frontal white matter.  Given our
% anatomical/ROI resources with the DK2009 atlas, we have to get a bit
% creative

% first we take the superior frontal gyrus.  However, this roi extends
% down the anterior frontal region of the brain, which mean's its most
% inferior border will end up somewhat low in the brain.  This is
% unfortunate though, be cause the posterior portion of this ROI
% terminates at about the height we would want to limit our search at
% (given that the SLF 1, for example, would be traveling within the
% white matter of this gyrus.
    frontalGyrusROI=bsc_roiFromAtlasNums(atlas,116+sidenum,1);
    
% We'll use the anterior border of the thalamus to clip the frontal gyrus
% ROI
    antThalPlane=bsc_planeFromROI_v2(thalLut(leftright),'anterior',atlas);
    posteriorFrontalGyrus=bsc_modifyROI_v2(atlas,frontalGyrusROI,antThalPlane,'posterior');
% Now we use the clipped frontal gyrus ROI to establish a ceiling limit for
% the cingulum
    inferiorBorderPosteriorFrontalGyrus=bsc_planeFromROI_v2(posteriorFrontalGyrus,'inferior',atlas);
    

% we can reuse the posterior border we generated with thirdVentPosterior
% to clip the posterior border of this ceiling roi, as the cingulum does
% have fanning in the posterior.  We'll essentially be generating an L
% shaped ROI, mirrored and rotated 90% clockwise, to prevent streamlines
% from arching over our cut ceiling exclusion into the superior white matter.
% We first have to create its two planar components though.
    coronalExclusion=bsc_modifyROI_v2(atlas,thirdVentPosterior,inferiorBorderPosteriorFrontalGyrus,'superior');
    transverseExclusion=bsc_modifyROI_v2(atlas,inferiorBorderPosteriorFrontalGyrus,thirdVentPosterior,'anterior');

% Merge these two planar rois to create the exclusion plane
    superiorToCingExclusionROI=bsc_mergeROIs(coronalExclusion,transverseExclusion);
% now apply the exclusion
    [~, noSuperiortraversalBool] = wma_SegmentFascicleFromConnectome(wbfg, {superiorToCingExclusionROI}, {'not'}, 'arbitraryName');
    %--------------------------------------------------------------------------
    
    occipitalBorder=bsc_planeFromROI_v2(110+sidenum,'posterior',atlas);
    
    cingPosteriorMidpointsBool=bsc_applyMidpointCriteria(wbfg,occipitalBorder,'anterior',lateralThalPlane,'medial');
    
    
    cingAnteriorMidpoints=bsc_applyMidpointCriteria(wbfg,thirdVentPosterior,'anterior');
    
    palTopPlane=bsc_planeFromROI_v2(palLut(leftright),'superior',atlas);
    
      palTopMidpoints=bsc_applyMidpointCriteria(wbfg,palTopPlane,'superior');
    %--------------------------------------------------------------------------

    
    
    %lentiLut(leftright,:)
    frontalSubocrticalROI=bsc_roiFromAtlasNums(atlas,[nACLut(leftright) 132+sidenum],1);
    
    [~, frontalSubcortBool] = wma_SegmentFascicleFromConnectome(wbfg, {frontalSubocrticalROI}, {'endpoints'}, 'arbitraryName');
    

    %--------------------------------------------------------------------------
        [lateralThalPlane] =bsc_planeFromROI_v2(thalLut(leftright), 'lateral',atlas);
    

   eitherLatThalBool=bsc_applyEndpointCriteria(wbfg,lateralThalPlane,'medial','both');
   
   lateralMidpointBorder=bsc_planeFromROI_v2(106+sidenum, 'lateral',atlas);
    %--------------------------------------------------------------------------
    midlineborder=bsc_planeFromROI_v2(ventricleLut(leftright), 'medial',atlas);
    boundedMidpointsBool=bsc_applyMidpointCriteria(wbfg,midlineborder,'lateral',lateralMidpointBorder,'medial');
    
  
     %-------------------------------------------------------------------------- 
     
     
    tractNameVar=strcat(sideLabel{leftright},'_parieto-frontal_Cingulum');
    classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontalConjunctBool,eitherLatThalBool,cingAnteriorMidpoints,noSuperiortraversalBool,posteriorEndpointBool,anteriorEndpointBool,boundedMidpointsBool);
    %bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)
    
    tractNameVar=strcat(sideLabel{leftright},'_subcallosal_Cingulum');
    classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,subcallosalConjunctBool,palTopMidpoints,frontalSubcortBool,noSuperiortraversalBool,boundedMidpointsBool);
    %bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)

    
    %initialize classification structure
classificationOut=[];
classificationOut.names=[];
%make sure it is the same length as the input
classificationOut.index=zeros(length(wbfg.fibers),1);
    
    bsc_quickPlotClassByName(wbfg,categoryClassification,'pericallosal')
    

relevanIndexes=find(contains(categoryClassification.names,'hipp'));

bsc_quickPlotClassByName(wbfg,categoryClassification,relevanIndexes)
    
    
    
    
    
end