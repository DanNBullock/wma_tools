function bsc_classifiedStreamEndpointCortex_BL

%% loading from config
config = loadjson('config.json');

wbfg=wma_loadTck(config.track);

if isfield(config,'output')
    load(config.output)
    classification=classification;
elseif isfield(config,'classification')
    load(config.classification)
    classification=classification;
end

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
