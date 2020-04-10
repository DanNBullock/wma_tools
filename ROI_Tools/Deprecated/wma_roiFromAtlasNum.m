function mergedROI = wma_roiFromAtlasNum(pathToAtlas,ROInums,smoothKernel)
%mergedROI = wma_roiFromFSnums(fsDir,fsROInums, smoothFlag, smoothKernel)
%
% Build a VISTASOFT ROI from a FreeSurfer segmentations and a series of
% indeces to regions in the FreeSurfer segmentaion
%
%  INPUTS: -fsDir:  path to THIS SUBJECT'S Atlas, ideally in nii.gz format.
%  Will attempt to convert .mgz
%
%  -fsROInums:  a list of region ID numbers from atlas that you would like
%  merged into a single ROI.  i.e. [ 3, 4 , ... 41,42]
%
%  -smoothKernel: the magnitude of the desired smoothing kernel.  Default
%  if none is input is 0, which is different than the default for
%  wma_roiFromFSnums which is 3
%
%  OUTPUTS: -mergedROI: a merged ROI which has coordinates corresponding to
%  each of the atlas regions associated with a number entered in the
%  fsROInums.  For example if pointing to the appropriate Freesurfer Atlas
%  [3,42] would create a grey matter mask roi, while [2, 41] would create a
%  white matter mask roi
%
%  ***NOTE***:  Can make a call to mri_convert, so requires that FreeSurfer
%  be installed and set up properly.
%
% (C) Daniel Bullock 2018 Indiana University

%% obtain atlas File
if ~exist(pathToAtlas,'file')
    warning('%s not found, checking for .mgz',strcat(pathToAtlas)) 
    %apaprently necessary for matlab?
    spaceChar={' '};
    [atPath atName atFile]=fileparts(pathToAtlas)
    if exist(strcat(atPath,atName,'.mgz'),'file')
    cmndString=strcat('mri_convert',spaceChar,atPath,atName,'.mgz',spaceChar, atPath,atName,'.nii.gz');

    [status result] = system(cmndString{1},'-echo');
    if status~=0
        error('/n Error generating aseg nifti file.  There may be a problem finding the .mgz file.  Ensure mri_convert is loaded. Output: %s',result)
    end
    else
        error('/n Neither atlas nifti nor mgz file found.')  
    end
else 
   atlasNifti=niftiRead(pathToAtlas); 
end

%reads in label data


%% iteratively use aparc atlas to make combined ROI

% get size of atlasNifti.data and make a blank matrix mask for it
atlasDataSize = size(atlasNifti.data);
blankLabelNifti(1:atlasDataSize(1),1:atlasDataSize(2),1:atlasDataSize(3))=false;

ROImask=blankLabelNifti;
roiNameString=[];
% sequentially build/amalgamate ROI masks and corresponding name string
for iRois = ROInums
    ROImask=or(ROImask,atlasNifti.data==iRois);
    if iRois==ROInums(end)
        roiNameString=strcat(roiNameString,num2str(iRois));
    else
        roiNameString=strcat(roiNameString,num2str(iRois),'_');
    end
end

% converts nifti formatting to roi formatting
% maybe include as ouput option later

%smooth if you want to
if exist('smoothKernel')
   ROImask =~ ((smooth3(ROImask,'box', smoothKernel)) == 0 );
end

%get index coordinates from nifti mask volume
Roi = find(ROImask);
[x1,y1,z1] = ind2sub(size(atlasNifti.data), Roi);

%create new ROI structure
if exist('smoothKernel')
mergedROI = dtiNewRoi([atName, roiNameString,'smooth',num2str(smoothKernel)], 'r');
else
mergedROI = dtiNewRoi([atName, roiNameString], 'r');
%apply the qto_xyz transform to transform the mask coordinates to acpc
%cordinates
mergedROI.coords = mrAnatXformCoords(atlasNifti.qto_xyz, [x1,y1,z1]);
end

end