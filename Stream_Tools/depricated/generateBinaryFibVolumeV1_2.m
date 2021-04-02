function [outputNifti]=generateBinaryFibVolumeV1_2(voxelResize, refNifti, track, outputType ,threshold, smoothKernel)
% [fiberBoolNifti, fiberDensityNifti]=generateBinaryFibVolume(voxelResize, T1, fg, threshold, smoothBool)
%
% Takes in a collection of streamlines (as a fiber grou) and converts it to
% a resampled volume and then smooths the resultant volume object.  Outputs
% a nifti for both the voxel-wise density of the fibers and a boolean mask
% of the fiber volume.
%
%  INPUTS
%
%  voxelResize:  What size you would like the voxels of the output niftis
%  to be.  BE SURE TO ENTER THIS AS A TRIPPLE (e.g. [.5 .5 .5] to indicate
%  voxel size in all dimensions)
%
%  refNifti:  The reference nifti that you would like the output to be
%  aligned to.  It is assumed that the input fiber group and input nifti
%  are aligned to one another.
%
%  outputType:  a string varible indicating which of the following you
%  would like for the output%
%  
%    -- 'binary':  a binary mask (only ones and zeros) in the voxels
%  corresponding to the track, in accordance with the input parameters.
%  
%    -- 'node':  a nifti with voxel values corresponding to to the
%  number of nodes per voxel *normalized* by the maximum nodes per voxel
%  found in this
%  
%    -- 'stream':  a nifti with voxel values corresponding to to the
%  number of nodes per voxel *normalized* by the total number of
%  streamlines from the input 'tract' variable.
%
%  tack:  a track object (vistasoft format) that you would like to turn
%  into a volume
%
%  threshold:  either the minimum absolute number of nodes (for both 'node'
%  and 'stream' option) per voxel (below which will be dropped to zero) or
%  the proportion of the max density (below which will be dropped to zero).
%  If integer, will interpret as nodes per voxel, if some number between 0
%  and 1, will interpret this as a density threshold)
%
%  smoothKernel:  the *diameter* of the smoothing kernel to be used.  *is
%  associated with the voxelResize value chosen, NOT absolute mm*.  In
%  accordance with matlab smooth3 rules, must be an odd number.
%
%  OTUPUTS
%
%  outputNifti:  The output nifti in accordance to with the outputNiftiType
%  input variable ->
%
%    -- 'binary':  a binary mask (only ones and zeros) in the voxels
%  corresponding to the track, in accordance with the input parameters.
%  
%    -- 'node':  a nifti with voxel values corresponding to to the
%  number of nodes per voxel *normalized* by the maximum nodes per voxel
%  found in this
%  
%    -- 'stream':  a nifti with voxel values corresponding to to the
%  number of nodes per voxel *normalized* by the total number of
%  streamlines from the input 'tract' variable.
%  
%
% First Vers 1.2 (3/14/17)
%  Refactored 4/26/2020
%% begin code


%convert the track's nodes to a coord cloud
inputFibers=track.fibers;
%create blank coordinate vector
coordVec=[];
for iFibers=1:length(inputFibers)
    coordVec=horzcat(coordVec,inputFibers{iFibers});
end

%turn those coords into an roi
%which way is this supposed to be oriented?
trackROI = dtiNewRoi(track.name, 'r', coordVec');

[trackNifti, ~] = dtiRoiNiftiFromMat(trackROI,refNifti,track.name,false);

[resizeNifti, ~, ~] = bsc_ResampleNifti(trackNifti, voxelResize, [], []);

keyboard




fiberBoolNifti=T1;
imgDim=T1.dim(1:3);
imgRes=T1.pixdim(1,1);
imgResize=imgDim*imgRes;
imgBins=ceil(imgResize/voxelResize);
emptyMatrix=zeros(imgBins);


for ifibers=1:length(fg.fibers)
    
    roundedExpandedTract=round(fg.fibers{ifibers}*(imgRes/voxelResize));
    
    for inodes=1:length(fg.fibers{ifibers})
        emptyMatrix(  roundedExpandedTract(1,inodes),roundedExpandedTract(2,inodes),roundedExpandedTract(3,inodes))= ...
            emptyMatrix(  roundedExpandedTract(1,inodes),roundedExpandedTract(2,inodes),roundedExpandedTract(3,inodes)) +1;
    end
end

%adjusts nifti object field information (probably not necessary)
fiberBoolNifti.dim=imgBins;
fiberBoolNifti.pixdim=[voxelResize voxelResize voxelResize];
fiberBoolNifti.qoffset_x=fiberBoolNifti.qoffset_x*(imgRes/voxelResize);
fiberBoolNifti.qoffset_y=fiberBoolNifti.qoffset_y*(imgRes/voxelResize);
fiberBoolNifti.qoffset_z=fiberBoolNifti.qoffset_z*(imgRes/voxelResize);

%duplicates object
fiberDensityNifti=fiberBoolNifti;

%% compute DiffusivityRatio of fg
% The "diffusivity ratio" of a fg (not a technical term, has nothing to do
% with the FA or anything) is the ratio of voxels with more than one node
% in them to the total number of non zero voxels.  Presumably, this is a
% measure of how "dense" the fg is.  If there are very few voxels with
% multiple nodes then this means that there may not be a very dense
% representation of the fiber.  Basically, I am developing this ratio as a
% potential benchmark for the viability of the volumes/surfaces we are
% developing.  My sense is that there may be some ratio value that would
% reliably lead to better printing outcomes.

moreThanOneNode=length(find(emptyMatrix>1));
nonZeroNode=length(find(emptyMatrix>0));
DiffusivityRatio=moreThanOneNode/nonZeroNode;


%% actually do smoothing if necessary
% smooths if necessary.  In theory de-islands as well, though I've never
% seen it in action.  It is set to go to keyboard if it detects islands.
% It could be that the function to detect islands is broken, so this may
% simply be pointless.
if smoothBool==0
    boolMatrixVersion=emptyMatrix>threshold;
    fiberBoolNifti.data=boolMatrixVersion;
else
    smoothData=smooth3(emptyMatrix,'gaussian',smoothKernel);
    % auto threshold computation
    % computes the appropriate threshold for the given percentile value.
    % Probably not computationally efficient, but principled (aside from the
    % arbitrary threshold) and adaptive to the specific case.
    uniqueVals=unique(smoothData);
    %calculate bin interval according to uniquevals
    binz=0:(uniqueVals(end)/10000):uniqueVals(end);
    %is there something wrong here with the calculation of the second bin?
    hisDist=histcounts(smoothData,binz);
    %count of nonzero entries, used for percentile calculation later
    nonZeroTotal=sum(hisDist(2:end));
    %loop to calculate cumulative sum of distribution
    for ibins=1:length(hisDist)
        if ibins==1
            cumulativeSums(ibins)=0;
        else
            cumulativeSums(ibins)= cumulativeSums(ibins-1)+hisDist(ibins);
        end
    end
    %calculate percentiles
    percentiles=(cumulativeSums/nonZeroTotal)*100;
    failureVec=(percentiles>thresholdPercent);
    %find corresponding bin value
    thresholdBinVal=binz(min(find(failureVec)));
    
    
    
    %the bwislands is currently returning an empty output for me, at least
    %on this test fiber.   I'm not sure if that means there are no islands
    %to be found for this fiber, or if the function doesn't work.  This
    %part may cause a problem at some point.
    islands=bwislands(smoothData);
    if ~isempty(islands)
        islandFlag=true;
        smoothData=deislands3d(smoothData,3);
    end
    
    boolMatrixVersion=smoothData>thresholdBinVal;
    fiberBoolNifti.data=boolMatrixVersion;
end
%convert to unit8 from bool (necessary for niftiSave function)
fiberBoolNifti.data=uint8(fiberBoolNifti.data);
%fiberdensityNifti may be useless at this point because it is unsmoothed,
%and the indexing done by boolmatrixversion is going to index a great deal
%of empty voxels in the sparser, unsmoothed emptyMatrix object.
%
%
%fiberDensityNifti.data=emptyMatrix(boolMatrixVersion);
%
% Update:  now just saving the smoothed object.
fiberDensityNifti.data=smoothData(boolMatrixVersion);
end
