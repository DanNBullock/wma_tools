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

%now fixed to work locally, just have config.json point to first csv in
%project directory.  ASSUMES BIDS/BRAINLIFE FORMAT
if ~isdeployed
    firstCsvPath=csvPaths{1};
    allDirs = split(firstCsvPath,'/');
    projDirInd=max(find(contains(allDirs,'proj')));
    dataDirInd=projDirInd+2;
    dataDirFull=char(allDirs{dataDirInd});
    endDataNameInd=strfind(dataDirFull,'id-');
    dataDirHeader=dataDirFull(1:endDataNameInd-2);
    slashFind=strfind(firstCsvPath,'/');
    %why?
    firstCsvPath=char(firstCsvPath);
    projPath=firstCsvPath(1:slashFind(projDirInd));
    dirContents=dir(projPath);
    subjects={dirContents(contains({dirContents(:).name},'sub')).name};
    for iSubjects=1:length(subjects)
        curSubjDir=fullfile(projPath,subjects{iSubjects});
        subDirContent=dir(curSubjDir);
        subjDataDir=subDirContent(contains({subDirContent(:).name},dataDirHeader)).name;
        csvPaths{iSubjects}=fullfile(curSubjDir,subjDataDir,'output_FiberStats.csv');
    end
else
    for iFiles= 1:length(csvPaths)
        subjects{iFiles}=config.x0x5F_inputs{iFiles}.meta.subject;
    end
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
