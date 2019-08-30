function bsc_SegROIfromPairStringList_BL()
% bsc_genNiftiROIfromStringList(feORwbfg,atlas,ROIstring, smoothKernel)
%
% Given a string list of rois (in the specified format) loops over
% the list and generates nifti ROI files for each 
%
%  INPUTS:
%  -feORwbfg: either a string or an object, to either a wbFG or an FE
%  structure
%
%  -atlas: path to an atlas or an atlas.
%
%  -ROIstring:  a string list of atlas based roi specifications from the
%  atlas that you would like merged into a single ROI.  eg: '2 56 30 54; 34
%  654 \n 25 45 56; 23 \n 456 34; 35 75'.  Must correspond to the values in
%  the atlas nifti itself.
%
%  OUTPUTS: 
% -classification: 
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
%  NOTE:  Makes a call to mri_convert, so requires that FreeSurfer be
%  installed and set up properly.
%
%  (C) Daniel Bullock 2018 Bloomington, Indiana
%% Begin code

 if ~isdeployed
    disp('adding paths');
     addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
     addpath(genpath('/N/u/brlife/git/jsonlab'))
     addpath(genpath('/N/u/brlife/git/vistasoft'))
      addpath(genpath('/N/u/brlife/git/wma_tools'))
          addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
 end

%config = loadjson('/N/dc2/projects/lifebid/HCP/Dan/GitStoreDir/ROIs2ROIsSegment/config.json');
config = loadjson('config.json');

feORwbfg = wma_loadTck(config.track);

ROIstring=config.roiPairs;
smoothKernel=config.smoothKernel;

if isfield(config,'atlas')
    atlas=config.atlas;
end

if isfield(config,'ROI')
      ROIdir=config.ROI;
end

%% gen ROI
fprintf('Generating ROIs for the following indicies: \n %s',ROIstring);
%just in case
ROIstring=strrep(ROIstring,'\n',newline);
ROIstring=strrep(ROIstring,';',newline);
stringCells = splitlines(ROIstring);
mkdir('roi/')

bsc_GenROIfromPairStringList_BL()

for iROIs=1:length(stringCells)
 mergedROIs{iROIs}=fullfile(pwd,'roi',strcat('ROI',num2str(iROIs),'.nii.gz'));
end
    [classification]=multiROIpairSeg(feORwbfg,mergedROIs);
    save('/classification/classification.mat','classification')
    fprintf('\n classification structure stored with %i streamlines identified across %i tracts',...
    sum(classification.index>0),length(classification.names))
    wma_formatForBrainLife_v2(classification,feORwbfg)
end
