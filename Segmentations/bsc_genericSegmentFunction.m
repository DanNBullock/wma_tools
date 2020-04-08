function [classificationOut] =bsc_genericSegmentFunction(wbfg, fsDir, categoryClassification)
% [classificationOut] = bsc_genericSegmentFunction(wbfg, fsDir, categoryClassification)
%
% This is the template segmentation script for the refactored version of
% wma tools segmentations.  This is being provided to serve as a basic
% framework for future segmentation functions and to illustrate the basic
% logic of this format.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to the freesurfer directory for the subject 
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

%set a path to the atlas you will be using.  In truth, you can use any
%atlas.  The key is that a number of functions (i.e. bsc_planeFromROI_v2 or
%bsc_roiFromAtlasNums
atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

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
    
    %% Tract 1
    %begin segmentation of tract 1
%=========================================================================   
    %1.  ESTABLISH CATEGORY CRITERIA
    %usingthe category segmentation output, go ahead and establish a
    %boolean index of the streamlines meeting this criteria.  Will be used
    %as part of a later logical conjunct.  Here we use the example of
    %within frontal lobe streamlines.
    frontoFrontalBool=  or( bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_frontal')),   bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_frontal_ufiber')));
%--------------------------------------------------------------------------  


    %2.  ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA
    %if you have more specific criteria (positive or negative) for the
    %cortical terminations of your tract, do that here.
%=========================================================================    
        %2.1 extract the relevant rois from the atlas
        %check function description for more details, leave final integer
        %input at 1 for no inflation
        [someROI1] =bsc_roiFromAtlasNums(atlasPath,[ROInums],1);
        [someROI2] =bsc_roiFromAtlasNums(atlasPath,[ROInums],1);
        
        %2.2  Use ROIs to find the indexes of streamlines which terminate
        %in those ROIS.  Consider also using bsc_endpointAtlasCriteria,
        %which does not require the rois to be extracted first, but does
        %take a while if you send in the whole brain fiber group.  It may
        %be possible to make a variant of that function which impliments a
        %speedup by preselecting using a bolean index input.
        [~, corticalCriteriaBool] =  bsc_tractByEndpointROIs(wbfg, {someROI1 someROI2});
    
        %2.3 Add aditional criteria if you wish, for example if you want 
        %endpoints to be above a particular anatomical roi.  Check the
        %documentation for these functions to see how to apply them with
        %other anatomical relations (e.g. 'medial', 'lateral', etc.)
        superiorPlaneROI = bsc_planeFromROI_v2([ROInums], 'superior',atlasPath);
        [relativeAnatomicalEndpointCriteriaBool]=bsc_applyEndpointCriteria(wbfg, superiorPlaneROI, 'superior','one');
%--------------------------------------------------------------------------


    %3.  APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
    %at this point in the segmentation, you've already established where
    %you want the streamlines to terminate, but this may still
    %underdetermine the tract of interest.  For example, I could be
    %interested in fronto occipital streamlines, but these could go via a
    %ventral (i.e. IFOF) or dorsal (i.e. the putative "SFOF", which probably
    %doesn't really exist).  To further subselect, we can apply additional
    %criteria.  We'll use an example of putting a plane above the posterior
    %limit  of the thalamus thalamus, such that we are (to some extent)
    %selecting for streamlines taking the 
%==========================================================================    
        %3.1  Application of anatomically informed planes
        [superiorThalPlane]= bsc_planeFromROI_v2(thalLut(leftright), 'superior',atlasPath);
        [posteriorThalPlane] = bsc_planeFromROI_v2(thalLut(leftright), 'posterior',atlasPath);
        %now that we have the two planes we use the function
        %bsc_modifyROI_v2 to use the superior plane to slice the posterior
        %plane.  Note: this would result in a different output if you
        %switched the input rois.  Check the function documentation for
        %more details.
        [superiorPosteriorThalPlane]=bsc_modifyROI_v2(atlasPath,posteriorThalPlane, superiorThalPlane, 'superior');
        %ALSO NOTE: this function is fairly versitile.  You could instead
        %input a roi number for the first roi input (input 2) and then cut
        %it using the second roi input (which could instead be a specific
        %3d coordinate rather than a full plane, if one wished).  In this
        %way you can further subsegment rois if they are too large for your
        %purposes.
        
        %now that you have the plane, you can use whatever roi criteria
        %function you prefer.  Here we use a modified version of
        %feSegmentFascicleFromConnectome.  This function variant likely
        %needs to be refactored and updated itself.  Also remember, this
        %could be applied as a exclusion criteria as well (via 'not')
        [~, superiorPosteriorThalCriteriaBool] = wma_SegmentFascicleFromConnectome(wbfg, {superiorPosteriorThalPlane}, {'and'}, 'arbitraryName');
%--------------------------------------------------------------------------        
        %3.2  Application of anatomically informed volumetric rois
        %One potential desired application might be that you wish to select
        %streamlines which travel along or near a particular gyrus.  We can
        %use wma tools to select these white matter volumes.
        
        %note that we are inflating the cortical rois so that they will
        %expand into the white matter.  The inflation kernel (last
        %variable) can be modified as needed.
        [inflatedROI1] =bsc_roiFromAtlasNums(atlasPath,[ROInums],5);
        [inflatedROI2] =bsc_roiFromAtlasNums(atlasPath,[ROInums],5);
        [wmROI] =bsc_roiFromAtlasNums(atlasPath,wmLut{leftright},1);
        
        %here we take the intersection of the rois
        [roi1WMintersection] = bsc_intersectROIs(inflatedROI1, wmROI);
        [roi2WMintersection] = bsc_intersectROIs(inflatedROI2, wmROI);
        
        %now we take the volumetric overlap of the wm areas
        [roi1and2WMintersection] = bsc_intersectROIs(roi1WMintersection, roi2WMintersection);
        
        %the resultant roi can now be used as a segmentation criteria
        [~, wmVolumeCriteriaBool] = wma_SegmentFascicleFromConnectome(wbfg, {roi1and2WMintersection}, {'and'}, 'arbitraryName');
%--------------------------------------------------------------------------        
        %3.3  Application of other anatomically informed criteria
        %One additional set of criteria that one might like to apply
        %relates to the midpoints of the streamlines for the tract.  For
        %example may want the midpoints to be superior (or medial, lateral,
        %etc.) to some plane.  We can do this as well.
        
        %say we wanted to require our midpoints to be above the thalamus.
        %We could use the planar roi generated in section 3.1 to do this
        [superiorThalMidpointCriteriaBool]=bsc_applyMidpointCriteria(wbfg, superiorThalPlane,'superior');       
%--------------------------------------------------------------------------         
    %4.  APPLY ALL BOOLEAN CRITERIA
        %now that we have obtained all of our desired criteria, in the form
        %of several boolean vectors, it's time to apply them as a conjunct
        %(many ands) to obtain a classification structure featuring this
        %tract
    tractNameVar=strcat(sideLabel{leftright},'TractName');
    classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,frontoFrontalBool,corticalCriteriaBool,relativeAnatomicalEndpointCriteriaBool,superiorPosteriorThalCriteriaBool,wmVolumeCriteriaBool,superiorThalMidpointCriteriaBool);

%% TRACT 2
% now you can apply a similar series of requirements as above in order to
% segment another tract within this function.  Why would you want to
% segment multiple tracts in one function?  In some cases you may be using
% the same roi multiple times in order to segment several tracts.  It would
% thus be useful to segment them in the same function and reuse the roi so
% that you don't have to keep regenerating it.  Also, it may be the case
% that once you segment a tract, you can then exclude those tracts from a
% subsequent segmentation.  In this way, you would be applying a criteria
% of [not a member of this previous tract] to some subsequent tract.

%SPECIAL NOTE:  Be sure to use bsc_concatClassificationCriteria after
%you've generated all of your boolean criteria for each tract.  If you
%don't do this, it won't be added to the classification structure that is
%eventually output from this function
end


end