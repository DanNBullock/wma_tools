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
    %begin segmentation of Uncinate Fasiculus
%=========================================================================   
    %1.  ESTABLISH CATEGORY CRITERIA
%--------------------------------------------------------------------------  
    %The unicnate is, very straightforwardly, a fronto temporal tract.
    frontoTemporalBool=  bsc_extractStreamIndByName(categoryClassification,strcat(sideLabel{leftright},'_frontal_to_temporal'));
%=========================================================================     

    %2.  ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA
    %in order to prevent certian suprrious streamlines (potentially
    %implausable, and certianly not part of the uncinate) we can include
    %some additional segmentation logic to exclude these streamlines *while
    %also not being so strict as to make the segmentation brittle*.  Two of
    %these endpoint  criteria will be necessary for the uncinate.  Given
    %that we require the posterior cluster of streamlines to be located in
    %the anterior temporal region we can implement criteria relative to the
    %dorso-ventral axis and the rostro-caudal axis
%=========================================================================  
     
%--------------------------------------------------------------------------  
    %dorso-ventral criteria:
    %in general, the endpoints of the uncinate are fairly low in the brain.
    %furthermore, we note that the arc of the uncinate occurs relatively
    %close to the amygdala.  Because the streamlines of the uncinate
    %approach from the ventral side of the amygdala (as they then proceed
    %anteriorly to the frontal lobes), and reach their apex at around the
    %top of the amygdala, it is necessarily the case that at least one
    %endpoint of these streamlines (namely the posterior/inferior
    %endpoint) is below the top of the amygdala.  Therefore, it is
    %logically necessary that it is *not* the case that both streamline
    %endpoints are superior to the top of the amygdala (or else it couldn't
    %engage in the arcing behavior).  Lets translate this into segmentation
    %logic.
    
    %begin by generating a plane from the top of the amygdala
    [amygdalaTopPlane] =bsc_planeFromROI_v2(amigLut(leftright), 'superior',atlas);
    
    %now we find all streamlines that have both streamline endpoints
    %above this plane (we'll negate this criteria later, to adhere to our
    %logic). Note, this is different than requiring neither streamline to
    %be above this plane.  This logical operation still permits maximally
    %one endpoint per streamline to be above this plane.
    bothAboveAmygBool=bsc_applyEndpointCriteria(wbfg,amygdalaTopPlane,'superior','both');
    
    %next we apply a similar bit of logic on the rostro-caudal axis
%--------------------------------------------------------------------------
    %rostro-caudal criteria:
    %in general, the endpoints of the uncinate are fairly *anterior* in the brain.
    %as we noted above, the arc of the uncinate occurs relatively
    %close to the amygdala. This is true of the rostro-caudal axis as well.
    %The streamlines of the uncinate have all begun to arc forward anterior
    %of the posterior of the border of the amygdala.  A consequence,
    %similar to the one noted above, is that it is *not* the case that both
    %endpoints are posterior of this posterior amygdala border.  Lets
    %translate this into segmentation logic.
    
    %begin by generating a plane from the posterior of the amygdala
    [amygdalaPosteriorPlane] =bsc_planeFromROI_v2(amigLut(leftright), 'posterior',atlas);
    
    %now we find all streamlines that have both streamline endpoints
    %posterior to this plane (we'll negate this criteria later, to adhere
    %to our logic). Note, this is different than requiring neither
    %streamline to be posterior to this plane.  This logical operation
    %still permits maximally one endpoint per streamline to be posterior to
    %this plane.
    bothPosteriorAmygBool=bsc_applyEndpointCriteria(wbfg,amygdalaPosteriorPlane,'posterior','both');
    
    %a third criteria specific to the midpoint and the rostro-caudal plane
    %will also be necessary.

%=========================================================================   

    %3.  APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
    %two additional (non-endpoint) anatomically informed criteria are
    %necessary.
    %the first of these is a reflection of the endpoint rostro-caudal
    %criteria, while the other is a reflection of the
    %category-segmentation.  Specifically, the category segmentation is
    %insufficient to isolate the uncinate.
    %Use 
    % bsc_quickPlotClassByName(wbfg,categoryClassification,'frontal_to_temporal')
    % to confirm this for yourself.  Note that the majority of the
    % non-uncinate fibers appear to be arcuate fibers .  As such we need to
    % apply additional criteria to exclude these arcuate-related
    % streamlines
%========================================================================= 
    %--------------------------------------------------------------------------
    %midpointrostro-caudal criteria:
    %if the arc of the uncinate occurs anterior to the posterior border of
    %the amygdala, it stands to reason that the midpoint of those
    %streamlines would also occur anterior to this border.  Even in the
    %case of more sigmoidally shaped streamlines (as opposed to u or j
    %shaped), which do occur and have their endpoints potentially posterior
    %to the posterior amygdala, the midpoint would still be anterior of
    %this border if the anterior endpoint is to terminate in the anterior
    %frontal lobes.  This can rather straightforwardly be translated into
    %segmentation logic.
    
    %we already have the amygdalaPosteriorPlane so we don't need to
    %generate it again

    %now we find all streamlines that the midpoint anterior of this plane.
    midpointAntOfPosteriorAmygBool=bsc_applyMidpointCriteria(wbfg,amygdalaPosteriorPlane,'anterior');
 %--------------------------------------------------------------------------   
  
    %subcortical structures like the thalamus (and amygdala) tend to be relatively
    %invariant in their location across subjects.  This is fortunate
    %for us, because it looks like all of the arcuate-like fibers
    %extend *past* the posterior of the thalamus.  Lets generate a
    %planar roi that we can use as an exclusion criterion.
    [posteriorThalPlane] =bsc_planeFromROI_v2(thalLut(leftright), 'posterior',atlas);
        
    %now lets use that plane as an exclusion criteria 
    [~, posteriorThalExcludedBool] = wma_SegmentFascicleFromConnectome(wbfg, {posteriorThalPlane}, {'not'}, 'arbitraryName');
 
%=========================================================================          
    %4.  APPLY ALL BOOLEAN CRITERIA
        %now that we have obtained all of our desired criteria, in the form
        %of several boolean vectors, it's time to apply them as a conjunct
        %(many ands) to obtain a classification structure featuring this
        %tract
    tractNameVar=strcat(sideLabel{leftright},'_uncinate');
    classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,posteriorThalExcludedBool,frontoTemporalBool,~bothPosteriorAmygBool,~bothAboveAmygBool,midpointAntOfPosteriorAmygBool);    
%% Arcuate
%as we noted when segmenting the uncinate, the fronto temporal category
%seems to be primarily composed of uncinate and arcuate streamlines.  Lets
%use what we genrated for the uncinate to get the arcuate
%========================================================================= 
    %1.  ESTABLISH CATEGORY CRITERIA
    %The arcuate is also, very straightforwardly, a fronto temporal tract.
    %we dont need to regenerate the frontoTemporalBool, we still have it.
%=========================================================================       
    %3.  APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
    %we'll skip to using anatomically informed criteria for now, and come
    %back to endpoint-related criteria if we need to.
    %
    %as before, the category segmentation is insufficient to isolate the
    %arcuate.
    %as before, you can use 
    % bsc_quickPlotClassByName(wbfg,categoryClassification,'frontal_to_temporal')
    % to confirm this for yourself.  Note that the majority of the
    % non-uncinate fibers appear to be arcuate fibers (and vice-versa.  All
    % we have to do is negate the boolean vector we got for the posterior
    % thalamus exclusion plane segmentation.  In this way we will obtain
    % the indexes of streamlines that *do* pass through that posterior
    % plane.
    posteriorThalExcludedBool=~posteriorThalExcludedBool;


end


end