function [classification] =bsc_segmentAslant(wbfg, fsDir,varargin)
%
%[classification] =bsc_segmentAslant(wbfg, fsDir)
%
% This function automatedly segments the middle longitudinal fasiculus
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% -classification:  A standardly constructed classification structure
% identifying the streamlines corresponding to the left and right Aslant

% (C) Daniel Bullock, 2018, Indiana University
%% Begin code
categoryPrior=varargin{1};

%initialize classification structure
classification=[];
classification.names=[];
classification.index=zeros(length(wbfg.fibers),1);

%establish side labels for later name use
sideLabel={'left','right'};

atlasPath=fullfile(fsDir,'mri/aparc.a2009s+aseg.nii.gz');

%iterates through left and right sides
for leftright= [1,2]
    
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme for Destrieux 2009. left = 1, right = 2.  This will cause an
    %issue with the parietal ROI as we will see
    sidenum=10000+leftright*1000;
    
    %superior frontal ROI:
    %generates the roi for the superior frontal termination of the aslant.
    %We inflate in order to be fairly generous with terminations
    [superiorROI] =bsc_roiFromAtlasNums(atlasPath,[116]+sidenum,7);
    superiorROI.name='superior';
    
    % lateral frontal roi
    %generates the roi for the lateral frontal termination of the aslant.
    %We inflate in order to be fairly generous with terminations.  The
    %subtraction of 10000 is due to a switch from using the 2009 Destrieux atlas
    %to the 2006 DK atlas.  This is because we want the
    %frontal_inf-Triangular_part ROI in order to be in keeping with the
    %litterature on the Aslant
    [lateralROI] =bsc_roiFromAtlasNums(atlasPath,[114 112]+sidenum,1);
    lateralROI.name='lateral';
    
    %the inflated lateralROI includes some areas near the insula that
    %shouldn't be included.  This border is used to demarcate those regions
    latBorder = bsc_planeFromROI_v2([104]+sidenum, 'medial',atlasPath);
    
    %here we modify the lateralROI to remove those parts that are medial to
    %this border
    [lateralROI]=bsc_modifyROI_v2(atlasPath,lateralROI, latBorder, 'lateral');
    
    %segment those streamlines with endpoints in the aforemenitoned ROIs
    [~, keep]=bsc_tractByEndpointROIs(wbfg,[{superiorROI}, {lateralROI}]);
    
    %create a not ROI to prevent crossing streamlines
    [notHemi]=bsc_makePlanarROI(atlasPath,0, 'x');
 
    
    %institute an anterior border based on the pericallosal sulcus
    antPeriCall = bsc_planeFromROI_v2([167]+sidenum, 'anterior',atlasPath);
    
    %institute a posterior border based on the paracentral gyrus/sulcus
    antParSup=bsc_planeFromROI_v2([103]+sidenum, 'anterior',atlasPath);
    
    preventAnteriorBool= bsc_applyMidpointCriteria(wbfg,antParSup,'posterior');
    
    
    %find those streamlines which run afoul of one of the three criteria.
    %This is, in essence a set of additional exclusion criteria for
    %streamlines which are included in the initial segmentation, but
    %nonetheless do not meet criteria for being included in the tract
    [~,boundedInd]=wma_SegmentFascicleFromConnectome(wbfg, [{antPeriCall},{antParSup},{notHemi}], {'not', 'not','not'}, 'exclusions');
    
    %use the positive and negative criteria to make a classification
    %structure
 
    frontoFrontalBool=  or( bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_frontal')),   bsc_extractStreamIndByName(categoryPrior,strcat(sideLabel{leftright},'frontal_to_frontal_ufiber')));
    
    [classification]=bsc_concatClassificationCriteria(classification,strcat(sideLabel{leftright},'Aslant'),keep,boundedInd,frontoFrontalBool,preventAnteriorBool);
end

end