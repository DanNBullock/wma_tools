function bsc_MBAplotTracts_BL()
%bsc_MBAplotTracts_BL()
%
%Brainlife wrapper for bsc_plotClassifiedStreamsAdaptive_v2.  See original function for
%more details.

if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/u/brlife/git/encode'))
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
    addpath(genpath('/N/u/brlife/git/mba'))
    addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
end

config = loadjson('config.json');

if isfield(config,'fe')
    feORwbfg=config.fe;
    fibNum=length(feORwbfg.fg.fibers);
else
    feORwbfg=wma_loadTck(config.track);
    fibNum=length(feORwbfg.fibers);
end

mkdir(fullfile(pwd,'images'))
saveDir=fullfile(pwd,'images');

t1=niftiRead(config.t1);

if isfield(config,'output')
    load(config.output)
elseif isfield(config,'classification')
    load(config.classification)
end

if isfield(config,'subSelect')
    subSelect=config.subSelect;
    subSelect=str2num(subSelect);   
else
    subSelect=1:length(classification.names);
end
%%
%code stolen from format for brainlife
classificationGrouped=wma_classificationStrucGrouping(classification);
neededColors=length(classificationGrouped.names);
smallCM = distinguishable_colors(neededColors,'k');

%find names and appropriate order for tracts
for iTracts=1:length(classification.names)
nameList{iTracts}=classification.names{iTracts};
end
cm=[];
cm=zeros(length(classificationGrouped.names),3)
%create a color vector with color pairings in the correct locations
for iGroups=1:length(classificationGrouped.names)
    curIndexes=bsc_extractStreamIndByName(classificationGrouped,classificationGrouped.names{iGroups});
    curNames={classification.names{unique(classification.index(curIndexes))}};
    for iNames=1:length(curNames)
        namePlace=strcmp(curNames{iNames},nameList);
      cm(namePlace,:)=smallCM(iGroups,:);
    end
end  

colors=cm(subSelect,:);
%%
views={'saggital','coronal','axial'};

for figView=1:3
    bsc_plotClassifiedStreamsAdaptive_v2(feORwbfg, classification ,t1, views{figView}, saveDir,subSelect,colors)
end

end