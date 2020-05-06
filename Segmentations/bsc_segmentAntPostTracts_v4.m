function [classificationOut] =bsc_segmentAntPostTracts_v4(wbfg,atlas,categoryClassification)
% [classificationOut] = bsc_segmentAntPostTracts_v4(wbfg, atlas, categoryClassification)
%
% This function segments several anterior-posterior oriented tracts.
% Because many of the planes and anatomical references are shared amongst
% these tracts, they are segmented here together.

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

if ischar(wbfg)
wbfg=fgRead(wbfg);
else
    %do nothing
end

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


%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
%% Uncinate
% DOI
% begin segmentation of Uncinate Fasiculus

%=========================================================================   
% 1. UNCINATE - ESTABLISH CATEGORY CRITERIA
% The unicnate is, very straightforwardly, a fronto temporal tract.
%--------------------------------------------------------------------------  
% we can extract the indexes for the relevant streamlines from the
% actegory segmentation
    frontoTemporalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_temporal'));
%=========================================================================     

% 2. UNCINATE - ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA
% in order to prevent certian suprrious streamlines (potentially
% implausable, and certianly not part of the uncinate) we can include
% some additional segmentation logic to exclude these streamlines *while
% also not being so strict as to make the segmentation brittle*.  Two of
% these endpoint  criteria will be necessary for the uncinate.  Given
% that we require the posterior cluster of streamlines to be located in
% the anterior temporal region we can implement criteria relative to the
% dorso-ventral axis and the rostro-caudal axis. 
%--------------------------------------------------------------------------  
% dorso-ventral criteria:
% In general, the endpoints of the uncinate are fairly low in the brain.
% Furthermore, we note that the arc of the uncinate occurs relatively
% close to the amygdala.  Because the streamlines of the uncinate
% approach from the ventral side of the amygdala (as they then proceed
% anteriorly to the frontal lobes), and reach their apex at around the
% top of the amygdala, it is necessarily the case that at least one
% endpoint of these streamlines (namely the posterior/inferior
% endpoint) is below the top of the amygdala.  Therefore, it is
% logically necessary that it is *not* the case that both streamline
% endpoints are superior to the top of the amygdala (or else it couldn't
% engage in the arcing behavior).  Lets translate this into segmentation
% logic.
    
% Begin by generating a plane from the top of the amygdala
    [amygdalaTopPlane] =bsc_planeFromROI_v2(amigLut(leftright), 'superior',atlas);
    
% Now we find all streamlines that have both streamline endpoints
% above this plane (we'll negate this criteria later, to adhere to our
% logic). Note, this is different than requiring neither streamline to
% be above this plane.  This logical operation still permits maximally
% one endpoint per streamline to be above this plane.
    bothAboveAmygBool=bsc_applyEndpointCriteria(wbfg,amygdalaTopPlane,'superior','both');
    
% Next, we apply a similar bit of logic on the rostro-caudal axis
%--------------------------------------------------------------------------
% rostro-caudal criteria:
% In general, the endpoints of the uncinate are fairly *anterior* in the brain.
% as we noted above, the arc of the uncinate occurs relatively
% close to the amygdala. This is true of the rostro-caudal axis as well.
% The streamlines of the uncinate have all begun to arc forward anterior
% of the posterior of the border of the amygdala.  A consequence,
% similar to the one noted above, is that it is *not* the case that both
% endpoints are posterior of this posterior amygdala border.  Lets
% translate this into segmentation logic.
    
% Begin by generating a plane from the posterior of the amygdala
    [amygdalaPosteriorPlane] =bsc_planeFromROI_v2(amigLut(leftright), 'posterior',atlas);
    
% Now we find all streamlines that have both streamline endpoints
% posterior to this plane (we'll negate this criteria later, to adhere
% to our logic). Note, this is different than requiring neither
% streamline to be anterior to this plane.  This logical operation
% still permits maximally one endpoint per streamline to be posterior to
% this plane.
    bothPosteriorAmygBool=bsc_applyEndpointCriteria(wbfg,amygdalaPosteriorPlane,'posterior','both');
    
% A third criteria specific to the midpoint and the rostro-caudal plane
% will also be necessary.

%=========================================================================   

%3. UNCINATE - APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
% Two additional (non-endpoint), anatomically informed criteria are
% necessary. The first of these is related to the endpoint rostro-caudal
% criteria, while the other is retlated to  the category-segmentation.

%--------------------------------------------------------------------------

% midpointrostro-caudal criteria:
% If the arc of the uncinate occurs anterior to the posterior border of
% the amygdala, it stands to reason that the midpoint of those
% streamlines would also occur anterior to this border.  Even in the
% case of more sigmoidally shaped streamlines (as opposed to u or j
% shaped), which do occur and have their endpoints potentially posterior
% to the posterior amygdala, the midpoint would still be anterior of
% this border if the anterior endpoint is to terminate in the anterior
% frontal lobes.  This can rather straightforwardly be translated into
% segmentation logic.
    
% We already have the amygdalaPosteriorPlane so we don't need to
% generate it again

% Now we find all streamlines that the midpoint anterior of this plane.
    midpointAntOfPosteriorAmygBool=bsc_applyMidpointCriteria(wbfg,amygdalaPosteriorPlane,'anterior');
%--------------------------------------------------------------------------   

% posterior non-traversal criteria:
% As it should have been clear by now the category segmentation is
% insufficient to isolate the uncinate.  Use 
% bsc_quickPlotClassByName(wbfg,categoryClassification,'left_frontal_to_temporal')
% to confirm this for yourself.  Note that the majority of the
% non-uncinate fibers appear to be arcuate fibers.  How can we exclude
% these fibers?  By applying a plane that selectively targets the
% arcuate streamlines.
 
% Subcortical structures, like the thalamus (and amygdala), tend to be
% relatively invariant in their location across subjects.  This is
% fortunate for us, because it looks like all of the arcuate-like fibers
% extend *past* the posterior of the thalamus.  Lets generate a planar roi
% that we can use as an exclusion criterion.
    [posteriorThalPlane] =bsc_planeFromROI_v2(thalLut(leftright), 'posterior',atlas);
        
% Now lets use that plane as an exclusion criteria 
    [~, posteriorThalExcludedBool] = wma_SegmentFascicleFromConnectome(wbfg, {posteriorThalPlane}, {'not'}, 'arbitraryName');
 
%=========================================================================          
%4.  UNCINATE - APPLY ALL BOOLEAN CRITERIA
% Now that we have obtained all of our desired criteria, in the form
% of several boolean vectors, it's time to apply them as a conjunct
% (many ands) to obtain a classification structure featuring this

% frontoTemporalBool - category criteria 
% bothAboveAmygBool - dorso-ventral criteria
% bothPosteriorAmygBool - rostro-caudal criteria
% midpointAntOfPosteriorAmygBool - midpointrostro-caudal criteria
% posteriorThalExcludedBool - posterior non-traversal criteria

    tractNameVar=strcat(sideLabel{leftright},'_uncinate');
% note: given that we are applying a logical conjunct, it doesn't
% matter what order the boolean vectors are entered in.
    classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,posteriorThalExcludedBool,frontoTemporalBool,~bothPosteriorAmygBool,~bothAboveAmygBool,midpointAntOfPosteriorAmygBool);    
% Run this commented line to visualize
%bsc_quickPlotClassByName(wbfg,classificationOut,strcat(sideLabel{leftright},'_uncinate'))
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Arcuate
% As we noted when segmenting the uncinate, the fronto temporal category
% seems to be primarily composed of uncinate and arcuate streamlines.  Lets
% use what we genreated for the uncinate to get the arcuate
%========================================================================= 

%1.  ARCUATE - ESTABLISH CATEGORY CRITERIA
%The arcuate is also, very straightforwardly, a fronto temporal tract.
%we dont need to regenerate the frontoTemporalBool, we still have it.
%--------------------------------------------------------------------------
%its in this variable
%frontoTemporalBool
%=========================================================================     

%2. ARCUATE - ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA
% Our specication of the arcuate being a fronto-temporal tract has
% limited the terminations of streamlines to the temporal lobes.
% Theoretically, we could be more specific about which cortical areas we
% expect terminations in, but we don't actually need to be.  This allows
% the tractography to have more influence in the outcome of the
% segmentation.  That being said, there may still be the need to exclude
% certian endpoint areas that are simply inconsistent with what is known
% about the arcuate.  We'll do that in this section
%--------------------------------------------------------------------------
% medial-exclusion criteria
% One thing you may also notice (if you plot the relevant classification
% category) is that we appear to be getting some endpoints that are quite
% close to the midline.  This isn't consitent with characterizations of the
% arcuate, but may occur with probabilistic tractography.
       
% To adress this lets take a reliable subcortical structure, the
% thalamus, which sits fairly centrally in the brain, and has a lateral
% border that could help us differentiate between the more lateral
% arcuate and the more medial cingulum
    [lateralThalPlane] =bsc_planeFromROI_v2(thalLut(leftright), 'lateral',atlas);
    
% Now we find all streamlines that have both streamline endpoints
% medial to this plane (we'll negate this criteria later, to adhere
% to our logic). This is what we would expect of the cingulum, but not
% the arcuate.  This border is quite medial, and so requiring that that
% neither endpoint be medial to this is not particularly stringent and
% is entirely consitent with existing understanding of the arcuate.
   eitherLatThalBool=bsc_applyEndpointCriteria(wbfg,lateralThalPlane,'medial','both');
%--------------------------------------------------------------------------   
% posterior termination exclusion criteria:
% In a similar sense, we want to get rid of streamines which may still be
% fronto-temporal (and meet our other criteria) but none the less still
% end far too posterior/rostrally to be plausable candidates for the
% arcuate.  We can use the posterior border of the venticle as our
% landmark.  One might consider using the posterior border of the lateral
% fisure ([11/12]141 in this atlas), however, given that the arcuate arcs
% around this (and thus its body necessarily extends posterior to it)
% this might, in practice, be a somewhat stringent criteria.  Instead, we
% use the posterior border of the lateral ventricle given that this is
% far more posterior, and thus shouldn't be particularly stringent.  
   [postVentPlane] =bsc_planeFromROI_v2(ventricleLut(leftright), 'posterior',atlas);
   
% now find any streamlines which have endpoints posterior to this.  We'll
% negate this later.  This would be equivalent to requiring 'both' to be
% 'anterior'.
   eitherPostVentBool=bsc_applyEndpointCriteria(wbfg,postVentPlane,'posterior','either');

%=========================================================================  
%3. ARCUATE - APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
% Given the general arcing morphology of the arcuate, there are some
% general, non endpoint-related criteria that we can make.
%-------------------------------------------------------------------------- 
% posterior traversal criteria:
% As before, the category segmentation is insufficient to isolate the
% arcuate. also as before, you can use:
% bsc_quickPlotClassByName(wbfg,categoryClassification,'frontal_to_temporal')
% to confirm this for yourself.  Note that the majority of the
% non-uncinate fibers appear to be arcuate fibers (and vice-versa).  All
% we have to do is negate the boolean vector we got for the posterior
% thalamus exclusion plane segmentation.  In this way we will obtain
% the indexes of streamlines that *do* pass through that posterior
% plane.
    
%lets do that now
    posteriorThalIncludeBool=~posteriorThalExcludedBool;
    
%-------------------------------------------------------------------------- 
% inferior traversal exclusion criteria:   
% In certian cases, with probabalistic tractography it is possible to
% obtain (a slight few) streamlines which follow a path similar to the
% IFOF (or reflect the diffusion algorithm being influenced by this
% adjacent structure) but nonetheless terminate in the frontal and
% temporal lobes.  To eliminate these we can apply a midpoint criteria,
% and require that our streamlines have midpoints above the top of the
% paladium, which would be relatively superior to the midpoints of any
% IFOF-like streamlines, but well below the main body or arch of the arcuate.
   [palTopPlane] =bsc_planeFromROI_v2(palLut(leftright), 'superior',atlas);
   midpointSupOfPalBool=bsc_applyMidpointCriteria(wbfg,palTopPlane,'superior');
   
%=========================================================================  
%4.  ARCUATE - APPLY ALL BOOLEAN CRITERIA
% Now that we have obtained all of our desired criteria, in the form
% of several boolean vectors, it's time to apply them as a conjunct
% (many ands) to obtain a classification structure featuring this  

% frontoTemporalBool - category criteria 
% eitherLatThalBool - medial-exclusion criteria
% eitherPostVentBool - posterior termination exclusion criteria
% posteriorThalIncludeBool - posterior traversal criteria
% midpointSupOfPalBool - inferior traversal exclusion criteria

   tractNameVar=strcat(sideLabel{leftright},'_arcuate');
% As before, order doesn't matter
   classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,~posteriorThalExcludedBool,frontoTemporalBool,~eitherLatThalBool,~eitherPostVentBool,midpointSupOfPalBool);
% Run this commented line to visualize
%bsc_quickPlotClassByName(wbfg,classificationOut,strcat(sideLabel{leftright},'_arcuate'))
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
%========================================================================= 
%%  IFOF
%  The IFOF is a anterior-posterior tract that runs inferior to the
%  arcuate, but has terminations in the frontal lobe (just as the uncinate
%  and arcuate do
%========================================================================= 
%1.  IFOF - ESTABLISH CATEGORY CRITERIA
% The IFOF is the paradigmatic (and perhaps only) example of a
% fronto-occipital track
%--------------------------------------------------------------------------

frontoOccipitalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_occipital'));
%=========================================================================  

%2.  IFOF - ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA
% Given the lack of other major strufrontoParietalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_parietal'));

   frontalGyrusROI=bsc_roiFromAtlasNums(atlas,116+sidenum,1);
[~, frontalGyrusBool] = wma_SegmentFascicleFromConnectome(wbfg, {frontalGyrusROI}, {'endpoints'}, 'arbitraryName');

antThalPlane=bsc_planeFromROI_v2(thalLut(leftright),'anterior',atlas)
posteriorFrontalGyrus=bsc_modifyROI_v2(atlas,frontalGyrusROI,antThalPlane,'posterior');

inferiorBorderPosteriorFrontalGyrus=bsc_planeFromROI_v2(posteriorFrontalGyrus,'inferior',atlas)

thirdVentPosterior=bsc_planeFromROI_v2(14,'posterior',atlas)

coronalExclusion=bsc_modifyROI_v2(atlas,thirdVentPosterior,inferiorBorderPosteriorFrontalGyrus,'superior')
transverseExclusion=bsc_modifyROI_v2(atlas,inferiorBorderPosteriorFrontalGyrus,thirdVentPosterior,'anterior')
superiorToCingExclusionROI=bsc_mergeROIs(coronalExclusion,transverseExclusion)

[~, noSuperiortraversalBool] = wma_SegmentFascicleFromConnectome(wbfg, {superiorToCingExclusionROI}, {'not'}, 'arbitraryName');

cingAnteriorMidpoints=bsc_applyMidpointCriteria(wbfg,thirdVentPosterior,'anterior')

periCTop=bsc_planeFromROI_v2(167+sidenum,'superior',atlas);
periClatBorder=bsc_planeFromROI_v2(167+sidenum,'lateral',atlas)
bsc_modifyROI_v2

lowMidpointsBool=bsc_applyMidpointCriteria(wbfg,periCTop,'inferior');

tractNameVar=strcat(sideLabel{leftright},'_SLF_Lat');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoParietalBool,bothLatThalBool);
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)

tractNameVar=strcat(sideLabel{leftright},'_SLF_Med');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoParietalBool,eitherLatThalBool,frontalGyrusBool,noSuperiortraversalBool);

tractNameVar=strcat(sideLabel{leftright},'_SLF_Cing');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoParietalBool,eitherLatThalBool,~cingAnteriorMidpoints,noSuperiortraversalBool);
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)ctures which exibit similar
% morphologies or trajectories to the IFOF, we do not need to specify any
% more specific endpoint criteria
%=========================================================================  

%3. IFOF - APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
% Here, we're going to use a somehat more advanced strategy.  One of the
% most characteristic features of the IFOF is the 'dip' that occurs
% adjacent to the lenticular nucleus.  Overall, our goal will be to apply
% an exclusion plane to prevent streamlines from traversing a particular
% region of the brain.  In order to enshrine that as a criteria we have
% approach it in several steps.
%--------------------------------------------------------------------------
% lenticular-dip criteria:
% First we have to specify the anterior-posterior location of where we want
% to do this.  The dip is near to its lowest takes place around the middle
% (in the anterior-posterior axis) of the lenticular nucleus, which, as it
% turns out, is about the anterior border of the thalamus. .
posteriorLentiPlane=bsc_planeFromROI_v2(thalLut(leftright),'anterior',atlas);

% We can reuse the plane generated earlier for the top of the globus paladus to
% specifiy the height we'd like to apply this plane at
% palTopPlane

% Now we cut the anterior-posterior plane (posteriorLentiPlane), such that we only take the part
% that is superior to the top of the globus paladus (palTopPlane)
superiorExclusionCutPlane=bsc_modifyROI_v2(atlas,posteriorLentiPlane, palTopPlane, 'superior');

% Now we can apply the exclusion plane
[~, superiorExclusionBool] = wma_SegmentFascicleFromConnectome(wbfg, {superiorExclusionCutPlane}, {'not'}, 'arbitraryName');
%--------------------------------------------------------------------------
% ventral traversal criteria:
% to invert the logic we applied to the arcuate, we can negate the
% midpointSupOfPalBool variable to require our streamlines to have a
% midpoint that is lower in the brain (as is appropriate for the IFOF)
% midpointSupOfPalBool

%=========================================================================  
%4.  IFOF - APPLY ALL BOOLEAN CRITERIA
% Now that we have obtained all of our desired criteria, in the form
% of several boolean vectors, it's time to apply them as a conjunct
% (many ands) to obtain a classification structure featuring this  
%--------------------------------------------------------------------------
% frontoOccipitalBool - category criteria 
% superiorExclusionBool - lenticular-dip criteria
% midpointSupOfPalBool - ventral traversal criteria

   tractNameVar=strcat(sideLabel{leftright},'_IFOF');
   classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoOccipitalBool,superiorExclusionBool,~midpointSupOfPalBool);
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
%=========================================================================  
%% SFOF
% The existance of a corresponding superior occipital fasiculus has been a
% subject of debate.  Regardless of the actual anatomical validity of this
% structure, we can nonetheless segment those streamlines (should any exist
% in an input connectome) that would constitute this putative tract.  Our
% abilitiy to segment a tract from a whole brain connectome is not to be
% taken as evidence for the existance of a tract.  Creedance for the
% existance of a tract can only be conferred relative to the degree to
% which the source whole brain connectome accurately recreates the
% connective archetecture of the brain.
%========================================================================= 
%1.  SFOF - ESTABLISH CATEGORY CRITERIA
% The SFOF, like the IFOF, is a fonto occipital tract (hence the name)
%--------------------------------------------------------------------------
% category criteria
% We don't need to reacquire this category boolean vector, as we already
% have it from the IFOF segmentation
%frontoOccipitalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_occipital'));
%=========================================================================  
%2.  SFOF - ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA 
%--------------------------------------------------------------------------
% medial-lateral division criteria
% In order to separate the more coherent lateral aspect of the putative
% SFOF, which would traverse a similar path as parts of the SLF x, from more
% medial streamlines, which would traverse a similar path to SLF x or the
% cingulum, we can use the lateral border of the thalamus once more
 bothLatThalBool=bsc_applyEndpointCriteria(wbfg,lateralThalPlane,'lateral','both');
%=========================================================================  
%4.  SFOF - APPLY ALL BOOLEAN CRITERIA
% Now that we have obtained all of our desired criteria, in the form
% of several boolean vectors, it's time to apply them as a conjunct
% (many ands) to obtain a classification structure featuring this  
%--------------------------------------------------------------------------

% frontoOccipitalBool - category criteria
% bothLatThalBool - medial-lateral division criteria

tractNameVar=strcat(sideLabel{leftright},'_Putative_Lateral_SFOF');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoOccipitalBool,midpointSupOfPalBool,bothLatThalBool);

tractNameVar=strcat(sideLabel{leftright},'_Putative_Medial_SFOF');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoOccipitalBool,midpointSupOfPalBool,~bothLatThalBool);
 %bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)
 %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 %% ILF
 
%=========================================================================  
%1.  ILF - ESTABLISH CATEGORY CRITERIA
%--------------------------------------------------------------------------
% category criteria
% The ILF is straightforwardly a occiptio temporal tact

occipitoTemporalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_occipital_to_temporal'));
%=========================================================================  
%2.  ILF - ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA 
% In order to separate the more coherent lateral aspect of the putative
% SFOF, which would traverse a similar path as parts of the SLF x, from more
% medial streamlines, which would traverse a similar path to SLF x or the
% cingulum, we can use the lateral border of the thalamus once more
%--------------------------------------------------------------------------
% posterior termination criteria
% The post central cingular gyrus tends to be a good indicator of the
% furtherst extent of the middle (as opposed to the further reaching
% inferior sections) occipital lobe.  This will be as far forward as we
% want our posterior endpoints proceeding
posteriorCingGyPlane=bsc_planeFromROI_v2(110+sidenum,'posterior',atlas);
% Although it may not be clear how, we will use this as a negative
% criteria.  We don't want streamlines that are extermely short, connecting
% to the medial temporal, lingual and fusiform gyri (such connections
% aren't consistent with traditional reports of the ILF).  Given that we
% are already requiring the streamlines to be Occipito-temporal, we can
% exclude streamlines that are too far anterior, because they would have
% both endpoints anterior of the plane we just made.
bothAntBool=bsc_applyEndpointCriteria(wbfg,posteriorCingGyPlane,'anterior','both');
%--------------------------------------------------------------------------
% anterior termination criteria
% To further be consistent with established reports, we'll require one
% endpoint to be sufficiently anterior to extend past the posterior border
% of the midbrain (fs region DC)

posteriorMBPlane=bsc_planeFromROI_v2(DCLut(leftright),'posterior',atlas);
antDCBool=bsc_applyEndpointCriteria(wbfg,posteriorMBPlane,'anterior','one');
%=========================================================================  
%3. ILF - APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
% In some tractography samples we get spurrious fibers traveling underneath
% the thalamus.  These are not consistent with established reports (as they
% are far too medial) and so we we will use a plane based on the thalamus
% to exclude these.
%--------------------------------------------------------------------------
% inferior-medial exclusion criteria
% We can use the posterior thalamus plane we made earlier for the uncinate
% and the lateral thalamus plane we made for the arcuate:
% posteriorThalPlane lateralThalPlane

% And clip this at the medial border of the thalamus
subThalPreventPlane=bsc_modifyROI_v2(atlas, posteriorThalPlane, lateralThalPlane, 'medial');

[~, inferiorExclusionBool] = wma_SegmentFascicleFromConnectome(wbfg, {subThalPreventPlane}, {'not'}, 'arbitraryName');
%=========================================================================  
%4.  ILF - APPLY ALL BOOLEAN CRITERIA
% Now that we have obtained all of our desired criteria, in the form
% of several boolean vectors, it's time to apply them as a conjunct
% (many ands) to obtain a classification structure featuring this 

% occipitoTemporalBool - category criteria
% bothAntBool - posterior termination criteria
% antDCBool - anterior termination criteria
% inferiorExclusionBool - inferior-medial exclusion criteria

tractNameVar=strcat(sideLabel{leftright},'_ILF');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,occipitoTemporalBool,~bothAntBool,antDCBool,inferiorExclusionBool);
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)
 %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
%%  SLF
   
%=========================================================================  
%1.  SLF - ESTABLISH CATEGORY CRITERIA
%--------------------------------------------------------------------------
% category criteria
% SLF 1 2 and 3 are all fronto-parietal tracts
frontoParietalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_parietal'));

   frontalGyrusROI=bsc_roiFromAtlasNums(atlas,116+sidenum,1);
[~, frontalGyrusBool] = wma_SegmentFascicleFromConnectome(wbfg, {frontalGyrusROI}, {'endpoints'}, 'arbitraryName');

posteriorFrontalGyrus=bsc_modifyROI_v2(atlas,frontalGyrusROI,antThalPlane,'posterior');

inferiorBorderPosteriorFrontalGyrus=bsc_planeFromROI_v2(posteriorFrontalGyrus,'inferior',atlas)

thirdVentPosterior=bsc_planeFromROI_v2(14,'posterior',atlas)

coronalExclusion=bsc_modifyROI_v2(atlas,thirdVentPosterior,inferiorBorderPosteriorFrontalGyrus,'superior')
transverseExclusion=bsc_modifyROI_v2(atlas,inferiorBorderPosteriorFrontalGyrus,thirdVentPosterior,'anterior')
superiorToCingExclusionROI=bsc_mergeROIs(coronalExclusion,transverseExclusion)

[~, noSuperiortraversalBool] = wma_SegmentFascicleFromConnectome(wbfg, {superiorToCingExclusionROI}, {'not'}, 'arbitraryName');

cingAnteriorMidpoints=bsc_applyMidpointCriteria(wbfg,thirdVentPosterior,'anterior')

periCTop=bsc_planeFromROI_v2(167+sidenum,'superior',atlas);
periClatBorder=bsc_planeFromROI_v2(167+sidenum,'lateral',atlas)
bsc_modifyROI_v2

lowMidpointsBool=bsc_applyMidpointCriteria(wbfg,periCTop,'inferior');

tractNameVar=strcat(sideLabel{leftright},'_SLF_Lat');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoParietalBool,bothLatThalBool);
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)

tractNameVar=strcat(sideLabel{leftright},'_SLF_Med');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoParietalBool,eitherLatThalBool,frontalGyrusBool,noSuperiortraversalBool);

tractNameVar=strcat(sideLabel{leftright},'_SLF_Cing');
classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoParietalBool,eitherLatThalBool,~cingAnteriorMidpoints,noSuperiortraversalBool);
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)
%bsc_quickPlotClassByName(wbfg,classificationOut,tractNameVar)

   %last criteria = sub lat hal exclude
   testFG=bsc_makeFGsFromClassification_v5(classificationOut,wbfg)
   bsc_plotROIandFG(testFG{1},periCTop,'r')
 
 
%initialize classification structure
classificationOut=[];
classificationOut.names=[];
%make sure it is the same length as the input
classificationOut.index=zeros(length(wbfg.fibers),1);




%% TODO

%aslant vert
%slf ant
%CC  cc
%cingulum  sub
% Thalamic  sub
% Parcvert 
% TPC  vert 
% MdLF -spl vert 
% MdLF -ang  vert
% VOF  vert
% OPTIC  sub
% Cerebellar  cere
% forceps Major  CC
% forceps Minor  CC

end


end