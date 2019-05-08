function wma_removeTractOutliers_BL()
% [classification] =wma_removeTractOutliers()
%
% This function removes tact streamline outliers in accordance with the user's input
% specifications

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% classificaton: a standardly organized classification structure
% (C) Daniel Bullock, Indiana University
%% Begin Code
if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
    addpath(genpath('/N/u/brlife/git/mba'))
    addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
end

config = loadjson('config.json');

if isfield(config,'output')
    load(config.output)
    classification=classification;
elseif isfield(config,'classification')
    load(config.classification)
    classification=classification;
end

centroidSD=config.centroidSD;
lengthSD=config.lengthSD;
maxIter=config.maxIter;

wbfg=wma_loadTck(config.track);

classification= removeOutliersClassification(classification,wbfg, centroidSD, lengthSD,maxIter);
mkdir(fullfile(pwd,'classification'));
save('classification/classification.mat','classification')

wma_formatForBrainLife_v2(classification,wbfg)
end




