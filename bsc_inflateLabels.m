function [inflatedAtlas] =bsc_inflateLabels(atlasPath,inflateItr)
%
% This function inflates the DK2009 atlas from freesurfer such that
% cortical labels are extended into the white matter and in to unknown
% labels
%
% Inputs:
% -inflateItr:  iterations of inflation.  Essentially akin to a radial
% smoothing kernel.
% -infTarg: label targets for inflation.  Default is wm, but can also
% include unlabeled and uknown.
% Outputs:
%  inflatedAtlas: the inflated version of the atlas.
% (C) Daniel Bullock, 2019, Indiana University

%atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

atlasIn=niftiRead(atlasPath);
fprintf('\n beginning island removal')
[ olab] = fnDeislandLabels_v2(atlasIn, [],5,999);

atlasIter=olab;
%set relevant ROI indicies
greyMatterROIS=[[101:1:175]+12000 [101:1:175]+11000];
subcorticalROIS=[10:13 17:20 26 58 27 49:56 59 ];
%spineROIS=[16 28 60];
%cerebellumROIS=[8 47 7 46 ];
%ventricleROIS=[31 63 11 50 4 43 77 14 24 15 44 5 62 30 80];
wmROIS=[41 2];
%ccROIS=[251:255];
%unknownROIS=[0 2000 1000];
%OpticCROI=[85];

%set number of neighboring voxels
threshold=6;
thresholdVal=threshold/27;

%set 
inflateTargROIs=[wmROIS 2000 1000 999];

%begin timer
tic

for iIter=1:inflateItr
    %initialize atlas objects
    atlasIterCur=atlasIter;
    atlasIterNext=atlasIter;
    
    %inflate the grey and subcortical rois
    inflateROis=bsc_roiFromAtlasNums(atlasIterCur,[greyMatterROIS subcorticalROIS], 3);
    
    %create nifti of the inflated rois
    [inflateInfo, ~] = dtiRoiNiftiFromMat(inflateROis,atlasPath,'roisInflate',false);
    
    %extract roi for target rois
    inflationTargets=bsc_roiFromAtlasNums(atlasIterCur,inflateTargROIs, 1);
    %report potential inflation target number
    fprintf('\n %i total voxels to inflate into',length(inflationTargets.coords))
    
    %create roi intersection of inflation targets and the inflated rois
    inflationMargin=bsc_intersectROIs(inflationTargets,inflateROis);
    %create a nii of the inflation marign, as a mask
    [ni, ~] = dtiRoiNiftiFromMat(inflationMargin,atlasPath,'atlasInflate',false);
    %smooth it, in order to find those voxels which are adjacent to the
    %mask
    niiSmoothedMask=smooth3(inflateInfo.data,'box',3)>=thresholdVal;
    %find the relevant indexes of the inflation margin(ni.data) which have
    %been inflated into (niiSmoothedMask)
    marginIndexes=find(ni.data&niiSmoothedMask);
    %get data size
    dataInSize=size(ni.data);
    %find indexes of this margin
    [marginX,marginY,marginZ]=ind2sub(dataInSize,marginIndexes);
    %determines if any of the inflating rois are on the edge of the nii
    XlimBool=or(marginX==1,marginX==dataInSize(1));
    YlimBool=or(marginY==1,marginY==dataInSize(2));
    ZlimBool=or(marginZ==1,marginZ==dataInSize(3));
    eliminationBool=any([XlimBool,YlimBool,ZlimBool],2);
    marginXMod=marginX(~eliminationBool);
    marginYMod=marginY(~eliminationBool);
    marginZMod=marginZ(~eliminationBool);
    
    fprintf('\n %i voxels for iteration %i',length(marginXMod),iIter)
    for iVoxels=1:length(marginXMod)
        
        candidateLabels=unique(atlasIterCur.data(marginXMod(iVoxels)-1:marginXMod(iVoxels)+1,marginYMod(iVoxels)-1:marginYMod(iVoxels)+1,marginZMod(iVoxels)-1:marginZMod(iVoxels)+1));
        candidateLabelsClean=setdiff(candidateLabels,inflateTargROIs);
        for iCandidates=1:length(candidateLabelsClean)
            labelCount(iCandidates)=sum(sum(sum(atlasIterCur.data(marginXMod(iVoxels)-1:marginXMod(iVoxels)+1,marginYMod(iVoxels)-1:marginYMod(iVoxels)+1,marginZMod(iVoxels)-1:marginZMod(iVoxels)+1)==candidateLabelsClean(iCandidates))));
        end
        
        tiebreakCount=1;
        if length(candidateLabelsClean(max(labelCount)==labelCount))>1
            
            %tiebreaker
            
            while length(candidateLabelsClean(max(labelCount)==labelCount))>1
                %will crap out of this happens on an edge
                candidateLabels=unique(atlasIterCur.data(marginXMod(iVoxels)-tiebreakCount:marginXMod(iVoxels)+tiebreakCount,marginYMod(iVoxels)-tiebreakCount:marginYMod(iVoxels)+tiebreakCount,marginZMod(iVoxels)-tiebreakCount:marginZMod(iVoxels)+tiebreakCount));
                candidateLabelsClean=setdiff(candidateLabels,inflateTargROIs);
                clear labelCount
                for iCandidates=1:length(candidateLabelsClean)
                    labelCount(iCandidates)=sum(sum(sum(atlasIterCur.data(marginXMod(iVoxels)-tiebreakCount:marginXMod(iVoxels)+tiebreakCount,marginYMod(iVoxels)-tiebreakCount:marginYMod(iVoxels)+tiebreakCount,marginZMod(iVoxels)-tiebreakCount:marginZMod(iVoxels)+tiebreakCount)==candidateLabelsClean(iCandidates))));
                end
                tiebreakCount=tiebreakCount+1;
            end
            
            %fprintf('\n tiebreaker resolved after %i iterations for %i voxels',tiebreakCount-1,length(find(atlasIterCur.data(voxelMask))))
            thisLabel=candidateLabelsClean(max(labelCount)==labelCount);
            % thisLabel=0;
            atlasIterNext.data(marginX(iVoxels),marginY(iVoxels),marginZ(iVoxels))=thisLabel;
            clear labelCount
        else
            thisLabel=candidateLabelsClean(max(labelCount)==labelCount);
            clear labelCount
            atlasIterNext.data(marginX(iVoxels),marginY(iVoxels),marginZ(iVoxels))=thisLabel;

            
        end
        
        
        
    end
    atlasIter=atlasIterNext;
end


ninesInd=find(atlasIter.data==999);
[xi yi zi]=ind2sub(size(atlasIter.data),ninesInd);

for Ireplace=1:length(xi)
    candidateLabels=unique(atlasIterCur.data(xi(Ireplace)-1:xi(Ireplace)+1,yi(Ireplace)-1:yi(Ireplace)+1,zi(Ireplace)-1:zi(Ireplace)+1));
    candidateLabelsClean=setdiff(candidateLabels,inflateTargROIs);
    if isempty(candidateLabelsClean)
        %artibrarily set it to the first one if there aren't any valid
        %targets
        candidateLabelsClean=candidateLabels(1);
    end
    
    for iCandidates=1:length(candidateLabelsClean)
        labelCount(iCandidates)=sum(sum(sum(atlasIterCur.data(xi(Ireplace)-1:xi(Ireplace)+1,yi(Ireplace)-1:yi(Ireplace)+1,zi(Ireplace)-1:zi(Ireplace)+1)==candidateLabelsClean(iCandidates))));
    end
    
    tiebreakCount=1;
    if length(candidateLabelsClean(max(labelCount)==labelCount))>1
        
        %tiebreaker
        
        while length(candidateLabelsClean(max(labelCount)==labelCount))>1
            %will crap out of this happens on an edge
            candidateLabels=unique(atlasIterCur.data(xi(Ireplace)-tiebreakCount:xi(Ireplace)+tiebreakCount,yi(Ireplace)-tiebreakCount:yi(Ireplace)+tiebreakCount,zi(Ireplace)-tiebreakCount:zi(Ireplace)+tiebreakCount));
            candidateLabelsClean=setdiff(candidateLabels,inflateTargROIs);
            clear labelCount
            for iCandidates=1:length(candidateLabelsClean)
                labelCount(iCandidates)=sum(sum(sum(atlasIterCur.data(xi(Ireplace)-tiebreakCount:xi(Ireplace)+tiebreakCount,yi(Ireplace)-tiebreakCount:yi(Ireplace)+tiebreakCount,zi(Ireplace)-tiebreakCount:zi(Ireplace)+tiebreakCount)==candidateLabelsClean(iCandidates))));
            end
            tiebreakCount=tiebreakCount+1;
        end
        
        %fprintf('\n tiebreaker resolved after %i iterations for %i voxels',tiebreakCount-1,length(find(atlasIterCur.data(voxelMask))))
        thisLabel=candidateLabelsClean(max(labelCount)==labelCount);
        % thisLabel=0;
        atlasIterNext.data(xi(Ireplace),yi(Ireplace),zi(Ireplace))=thisLabel;
        clear labelCount
    else
        thisLabel=candidateLabelsClean(max(labelCount)==labelCount);
        clear labelCount
        atlasIterNext.data(xi(Ireplace),yi(Ireplace),zi(Ireplace))=thisLabel;
        
        
    end
        atlasIter=atlasIterNext;
end
    
    
    fprintf('atlas inflation complete \n')
    toc
    
    inflatedAtlas=atlasIter;
    niftiWrite(inflatedAtlas,atlasPath)
    
end
