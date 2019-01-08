function bsc_feAndSegQualityCheck_BL()
%[figHandle, results]= bsc_feAndSegQualityCheck_BL(feORwbfg, classification, saveDir)
%
%Brainlife wrapper for bsc_feAndSegQualityCheck.  See original function for
%more details.


 if ~isdeployed
    disp('adding paths');
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

if isfield(config,'classification')
    load(config.classification)
      classification=classification;
end

saveDir=pwd;

bsc_feAndSegQualityCheck(feORwbfg, classification, saveDir)

end
