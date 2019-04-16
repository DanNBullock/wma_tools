function bsc_classifiedStreamEndpointCortex_BL
%bsc_classifiedStreamEndpointCortex(wbFG, classification, fsDir, saveDir,subSelect, decayFunc, decayRadiusThresh)

if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
    addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
end

%% loading from config
config = loadjson('config.json');

wbfg=wma_loadTck(config.track);

load(config.output)

if isfield(config,'subSelect')
    subSelect=config.subSelect;
else
    subSelect=1:length(classification.names);
end

saveDir=fullfile(pwd,'rois/rois/');
mkdir(saveDir);

fsDir=strcat(pwd,'/freesurfer');

decayFunc=config.decayFunc;

decayRadiusThresh=config.decayRadiusThresh;

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

if ~exist(atlasPath,'file')
    fprintf('\n FILE %i NOT FOUND',atlasPath)
end
%% endpointGeneration
bsc_classifiedStreamEndpointCortex(wbfg, classification, fsDir, saveDir,subSelect, decayFunc, decayRadiusThresh);

%% function