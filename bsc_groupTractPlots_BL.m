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


csvPaths=config.csv;

%will break if not on docker, due to strange reading behavior
 for iFiles= 1:length(csvPaths)
     %blcsvPaths{iFiles}=fullfile(fileDirs{iFiles},'output_FiberStats.csv');
     subjects{iFiles}=config.x0x5F_inputs{iFiles}.meta.subject;
 end

plotProperties=config.plotProperties;
if ~isempty(str2num(plotProperties))
    plotProperties=str2num(plotProperties);
else 
    plotProperties = erase(plotProperties,' ');
    plotProperties = strsplit(plotProperties,',');
    %do nothing, its presumably some strings
end
    
zThresh=config.zThresh;
zThresh=str2num(zThresh);

bsc_plotTractZscoreMeasures_pathsVersion(csvPaths,plotProperties,pwd)

if isfield(config,'zThresh')
bsc_saveTrackCheckList_pathsVersion(csvPaths,plotProperties,zThresh,subjects,pwd)
end

end
