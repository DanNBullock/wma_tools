function bsc_groupTractPlots_BL()
%% Begin Code
if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
    addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
end

config = loadjson('config.json');

fileDirs=config.csvDirs;

for iFiles= 1:length(fileDirs)
    csvPaths{iFiles}=fullfile(fileDirs{iFiles},'output_FiberStats.csv');
end

plotProperties=config.plotProperties;
if ~isempty(str2num(plotProperties))
    plotProperties=str2num(plotProperties)
else
    %do nothing, its presumably some strings
end
    
zThresh=config.zThresh;
zThresh=str2num(zThresh);

bsc_plotTractZscoreMeasures_pathsVersion(csvPaths,plotProperties,pwd)

if isfield(config,'zThresh')
bsc_saveTrackCheckList_pathsVersion(csvPaths,config.plotProperties,zThresh,pwd)
end

end
