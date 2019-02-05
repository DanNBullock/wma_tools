function bsc_feAndSegQualityCheck_BL()
%[figHandle, results]= bsc_feAndSegQualityCheck_BL(feORwbfg, classification, saveDir)
%
%Brainlife wrapper for bsc_feAndSegQualityCheck.  See original function for
%more details.


 if ~isdeployed
    disp('adding paths');
     addpath(genpath('/N/u/brlife/git/encode'))
     addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
     addpath(genpath('/N/u/brlife/git/jsonlab'))
     addpath(genpath('/N/u/brlife/git/vistasoft'))
     addpath(genpath('/N/u/brlife/git/wma_tools'))
 end

%config = loadjson('/N/dc2/projects/lifebid/HCP/Dan/GitStoreDir/ROIs2ROIsSegment/config.json');
config = loadjson('config.json');

if isfield(config,'fe')
    feORwbfg=config.fe;
else
    feORwbfg=config.wbfg;
end

saveDir=pwd;

if isfield(config,'output')
    load(config.output)
    classification=classification;
    bsc_feAndSegQualityCheck(feORwbfg, classification, saveDir)
else
    bsc_feAndSegQualityCheck(feORwbfg, [], saveDir)
end
end