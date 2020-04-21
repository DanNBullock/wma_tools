function [classificationOut] =bsc_segmentSuperficialFibers_v3(wbfg, atlas)
%[classificationOut] =bsc_segmentSuperficialFibers_v3(wbfg, atlas,inflateITer)
%
% This function automatedly segments out superficial fibers in a
% tractogram.  Here, we take supertifical to mean those streamlines which
% are close to the grey matter.  We use this (along with a lenght criteria)
% to select putative u-fibers
%
% Inputs:
% -wbfg: a whole brain fiber group structure
% -atlas: path to THIS SUBJECT'S %DK2009 Atlas in
% nii.gz format, or the object itself
%
%
% Outputs:
%  classificationOut:  standardly constructed classification structure,
%  indicating streamlines which are purported to be u-fibers
% (C) Daniel Bullock, 2019, Indiana University
% refactored 4/20/20

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

%we'll need to determine how deep into the white matter streamlines are
%going, so that means we'll need to be able to select wm voxels.  Lets go
%ahead and get those label indexes
wmLut=[2,41];

%for the same reason we'll also need the grey matter label indexes
greyMatterROIS=[[101:1:175]+11000; [101:1:175]+12000];

%we can go ahead and use this information to begin creating some initial
%ROIs.  Because some of these are whole brain, we can do them outside the
%left-right loop.  Before this though, lets inflate the cortical and
%subcortical labels of our brain, so that we can be generous with respect
%to streamlines that may be terminating near the GM/WM border.
%The terminal output of bsc_inflateLabels_v3 should tell you how many
%voxels were impacted by this inflation.
[inflatedAtlas] =bsc_inflateLabels_v3(atlas,2);
%we'll pass this atlas along in future steps

%lets use the grey matter labels to extract  left and right rois.
%Because our source for the GM roi is already inflated we don't need to
%inflate it again with our input parameters.
leftGMinf=bsc_roiFromAtlasNums(inflatedAtlas,greyMatterROIS(1,:),1);
rightGMinf=bsc_roiFromAtlasNums(inflatedAtlas,greyMatterROIS(2,:),1);


%one criteria that we may presume for u-fibers is that they run close to
%the edge of the white matter (adjacent to the grey matter, and thus quite
%superficial).  One method might be, for each streamline node, compute its
%distance from the grey matter border.  This would likely be
%computationally intensive.  An alternative huristic would be to require
%that the midpoint of each streamline be contained within that new inflated
%grey matter roi that we made (remember, we inflated it into the white
%matter)
%lets go ahead and merge our grey matter rois into a single roi that we can
%use to apply this heuristic
mergedGreyRoi=bsc_mergeROIs(leftGMinf,rightGMinf);
%and now lets apply the heuristic.
midpointBool=bsc_midpointROISegment(wbfg,mergedGreyRoi);
%now we have our first criteria for u-fibers.  For now, we're doing this
%outside the left-right loop because we wouldn't want to do it twice (by
%doing it within the loop)

%It would seem that another part of the definition of being a u-fiber is
%being relatively short.  Lets go ahead and compute all streamline lengths
%in the input wbfg
streamLengths = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), wbfg.fibers, 'UniformOutput', true);
%it's unclear what sort of length we should assume is an appropriate
%threshold for u-fibers.  For now we'll just arbitrarily select 30 mm
streamLengthThreshold=30;
%lets apply the threshold to the length measurements we just obtained
streamLengthsBool=streamLengths<streamLengthThreshold;
%now we have our second criteria

%for now, we've done all the work we can, outside of the left-right loop

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
    %% u-fibers
    %begin segmentation of u-fibers
%=========================================================================   
    %1.  ESTABLISH CATEGORY CRITERIA
%given that this is actually a prior for the category segmentation, we
%can't really presume a category for it now.  For now we'll leave this area
%blank
 %--------------------------------------------------------------------------  

    %2.  ESTABLISH MORE SPECIFIC ENDPOINT CRITERIA
    %minimally, we can say that u-fibers are corti-cortical, in that the
    %connect regions of the grey matter.  We can include this as a criteria
%=========================================================================    
        %lets apply a requirement that the streamlines for u-fibers have to
        %terminate (1) in the cortex and, (2) in the same hemisphere of the
        %cortex.  Note, this requirement for being non-interhemispheric is
        %a stipulation on our part.
        if leftright==1
        [~, corticalCriteriaBool] =  bsc_tractByEndpointROIs(wbfg, {leftGMinf leftGMinf});
        else
        [~, corticalCriteriaBool] =  bsc_tractByEndpointROIs(wbfg, {rightGMinf rightGMinf});
        end
        %this gives us our third criteria
%--------------------------------------------------------------------------


    %3.  APPLY GENERIC, ANATOMICALLY INFORMED CRITERIA
    %a futher anatomical prior we have for u-fibers is that they *do not*
    %travel in the deep white matter.  Lets impliment this as a
    %segmentation criteria
%==========================================================================
        %lets extract the WM roi from the inflated atlas.  Because of how
        %the grey matter has been inflated, the white matter we extract
        %from this atlas reprents deep white matter, away from the cortex.
        wmROIDeep=bsc_roiFromAtlasNums(inflatedAtlas,wmLut(leftright),1);
        
        %lets get those streamlines which *do not* intersect with this ROI,
        %in this way we are excluding all deep white matter streamlines &
        %tracts
        [~, superficialFibersBool]=wma_SegmentFascicleFromConnectome(wbfg, [{wmROIDeep}], {'not'}, 'dud');
        %Finally, this gives us our fourth criteria
         
%--------------------------------------------------------------------------         
    %4.  APPLY ALL BOOLEAN CRITERIA
        %now that we have obtained all of our desired criteria, in the form
        %of several boolean vectors, it's time to apply them as a conjunct
        %(many ands) to obtain a classification structure featuring this
        %tract
    tractNameVar=strcat(sideLabel{leftright},'_u-fibers');
    classificationOut=bsc_concatClassificationCriteria(classificationOut,tractNameVar,midpointBool,streamLengthsBool,corticalCriteriaBool,superficialFibersBool);


end %left-right loop complete

%segmentation complete
end