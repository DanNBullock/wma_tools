function densityMask = bsc_makeTractDensityMask(tract, referenceNifti, voxelResize, threshold, smoothParam, normalizeBool)
%densityMask = bsc_makeTractDensityMask(tract, voxelResize, threshold, smoothParam, normalizeBool)
%
% Purpose:  This function creates a density mask (NIfTI format) for an input
% tract.
%
% INPUTS
%
% tract: a tract structure containing some number of streamlines from which
% a streamline density mask will be derived
%
% referenceNifti:  The reference NIfTI that is to be used as the reference
% for generating the mask.  Theoretically it may be possible to write a
% function of this sort without using a reference NIfTI
%
% voxelResize:  the voxel resize parameter.  0 or [] results in no
% resizing.  Input value indicates the desired isometric voxel size of the
% output.  WARNING:  given that the node sampling rate is ~ 1mm, resizing
% of less than 1mm can result in masks with many holes.  Recomend use of
% smoothing to combat this.
%
% threshold:  The minimum number of streamlines (nodes, in practice) required in a voxel needed to
% prevent the value from being thresheld to 0.  This is applied AFTER the
% voxel resize but BEFORE any smoothing is applied.
%
% smoothParam: the smoothing kernel (necessarily odd, corresponding to
% diameter +1) to be applied to the mask.  Currently implemented as
% gaussian
%
% normalizeBool:  Whether or not to apply normalization.  This will occur
% at the last step.  Will divide by highest density such that values range
% from 0 to 1.
%
%  OUTPUTS
%
%  densityMask:  A NIfTI mask containing the density mask corresponding to
%  the input tract.
%
% Dan Bullock 2021
% Apparerently a repreated refactoring of the half complete generateBinaryFibVolumeV1_2
%
%% begin code

%extract the name
tractName=tract.name;

%perform resizing, if varible exists
if exist('voxelResize')
        if and(~isempty(voxelResize), ~voxelResize==0)
            [resizeNifti] = bsc_ResampleNifti(referenceNifti, [voxelResize voxelResize voxelResize], [], []);
        else
            resizeNifti=referenceNifti;
        end
else
    resizeNifti=referenceNifti;
end

%create new blank strcute for data
blankData=zeros(resizeNifti.dim);

fprintf('\n nifti resample complete, computing new streamline coorespondances')

%convert streams to image space
imgSpaceStreams=cellfun(@(x) mrAnatXformCoords(resizeNifti.qto_ijk,x), tract.fibers,'UniformOutput', false);

%memory saving technique
clear tract

%apply floor to streams, this gets us image space coords
roundedStreams=cellfun(@floor, imgSpaceStreams ,'UniformOutput', false);

%memory saving technique
clear imgSpaceStreams

%this eats all memory, dont do it.
%keeping it here as a proof of concept / indication of what is
%theoretically possible
%dataSum=cellfun(@(x) blankData(x(:,1),x(:,2),x(:,3))+1, roundedStreams ,'UniformOutput', false)

fprintf('\n streamline conversion complete, computing voxel traversal totals')

%can't think of elegant way to sum across streamlines
%loop across tract's streamlines
% for ifibers=1:length(roundedStreams)
%     ifibers
%     %extract current stream
%     currentFiber=roundedStreams{ifibers};
%     
%     %use the nodes to index into the blank mask and iterate by 1
%     blankData(currentFiber(:,1),currentFiber(:,2),currentFiber(:,3))=blankData(currentFiber(:,1),currentFiber(:,2),currentFiber(:,3))+1;
% end

%soichi's vastly faster version
for s = 1:length(roundedStreams)
    stream = roundedStreams{s};
    for p = 1:length(stream)
        x = stream(p,1);
        y = stream(p,2);
        z = stream(p,3);
        blankData(x, y, z) = blankData(x, y, z) + 1;
    end
end

fprintf('\n density masking operation complete')

%perform thresholding
if exist('threshold')
        if and(~isempty(threshold), ~threshold==0)
            blankData(blankData<threshold)=0;
        else
            %do nothing;
        end
else
    %do nothing;
end

%perform smoothing
if exist('smoothParam')
    if and(~isempty(smoothParam), ~smoothParam==0)
        %perform appropriateness check
        if isodd(smoothParam)
            blankData=smooth3(blankData,'gaussian',smoothParam);
        else
            error('input smoothParam is not odd, value must be odd to perform smoothing' )
        end
    else
        %do nothing;
    end
else
    %do nothing;
end

if exist('normalizeBool')
    if and(~isempty(normalizeBool), ~normalizeBool==0)
        blankData=blankData/max(max(max(blankData)));
    else
        %do nothing;
    end
else
    %do nothing;
end

fprintf('\n optional transformations complete.')

%set output NIfTI structure
densityMask=resizeNifti;

%set the data in the output NIfTI structure;
densityMask.data=blankData;
densityMask.fname=tractName;

fprintf('\n density mask creation for %s complete',tractName)

%function complete
end