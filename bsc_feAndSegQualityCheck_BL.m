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
    feORwbfg=dtiImportFibersMrtrix(config.track, .5);
end

saveDir=pwd;

if isfield(config,'output')
    load(config.output)
    classification=classification;
    bsc_feAndSegQualityCheck(feORwbfg, classification, saveDir)
else
    bsc_feAndSegQualityCheck(feORwbfg, [], saveDir)
end

load('tractomeResultStruc.mat')

tableArray{1,1}='wbfg';
tableArray{1,2}=results.WBFG.stream_count;
tableArray{1,3}=results.WBFG.volume;
tableArray{1,4}=results.WBFG.avg_stream_length;
tableArray{1,5}=results.WBFG.stream_length_stdev;
tableArray{1,6}=results.WBFG.avgFullDisp;
tableArray{1,7}=results.WBFG.stDevFullDisp;
tableArray{1,8}=results.WBFG.LogFitA;
tableArray{1,9}=results.WBFG.LogFitB;
tableArray{1,10}=results.WBFG.length_total;

%endpointDensity1
tableArray{1,11}=nan;
%endpointDensity2
tableArray{1,12}=nan;
%avgEndpointDist1
tableArray{1,13}=nan;
%avgEndpointDist2
tableArray{1,14}=nan;
%stDevEndpointDist1
tableArray{1,15}=nan;
%stDevEndpointDist2
tableArray{1,16}=nan;
%midpointDensity
tableArray{1,17}=nan;
%avgMidpointDist
tableArray{1,18}=nan;
%stDevMidpointDist
tableArray{1,19}=nan;
%norms.volumeProp
tableArray{1,20}=1;
%norms.countProp
tableArray{1,21}=1;
%norms.wireProp
tableArray{1,22}=1;

fieldNames={'stream_count', 'volume','avg_stream_length','stream_length_stdev','avgFullDisp','stDevFullDisp',...
    'LogFitA','LogFitB','length_total','endpointDensity1','endpointDensity2','avgEndpointDist1','avgEndpointDist2',...
    'avgEndpointDist2','stDevEndpointDist1','stDevEndpointDist2','endpointDensity1','endpointDensity2','midpointDensity',...
    'avgMidpointDist','stDevMidpointDist','norms.volumeProp','norms.countProp','norms.wireProp'};

fullFieldNames={'Tract Name','Streamline count', 'volume','avgerage streamline length','streamline length stdev','average full displacement','full displacement stdev',...
    'a*exp(b*x) fit a','a*exp(b*x) fit b','streamline length total','endpoint 1 density','endpoint 2 density','average endpoint distance from centroid 1',...
    'average endpoint distance from centroid 2','stdev of endpoint distance from centroid 1','stdev endpoint distance from centroid 2','midpoint density',...
    'average Midpoint distance from centroid','stDev of Midpoint distance from centroid','total volume proportion','total count proportion','total wiring proportion'};

conditionals=11:19;


if isfield(config,'output')
    for itracts=2:length(results.WBFG.tractStats)+1
        tableArray{itracts,1}=results.WBFG.tractStats{itracts}.name;
        tableArray{itracts,2}=results.WBFG.tractStats{itracts}.stream_count;
        tableArray{itracts,3}=results.WBFG.tractStats{itracts}.volume;
        tableArray{itracts,4}=results.WBFG.tractStats{itracts}.avg_stream_length;
        tableArray{itracts,5}=results.WBFG.tractStats{itracts}.stream_length_stdev;
        tableArray{itracts,6}=results.WBFG.tractStats{itracts}.avgFullDisp;
        tableArray{itracts,7}=results.WBFG.tractStats{itracts}.stDevFullDisp;
        tableArray{itracts,8}=results.WBFG.tractStats{itracts}.LogFitA;
        tableArray{itracts,9}=results.WBFG.tractStats{itracts}.LogFitB;
        tableArray{itracts,10}=results.WBFG.tractStats{itracts}.length_total;
        
        if ~isfield(results.WBFG.tractStats{itracts},'endpointDensity1')
            %endpointDensity1
            tableArray{itracts,11}=nan;
            %endpointDensity2
            tableArray{itracts,12}=nan;
            %avgEndpointDist1
            tableArray{itracts,13}=nan;
            %avgEndpointDist2
            tableArray{itracts,14}=nan;
            %stDevEndpointDist1
            tableArray{itracts,15}=nan;
            %stDevEndpointDist2
            tableArray{itracts,16}=nan;
            %midpointDensity
            tableArray{itracts,17}=nan;
            %avgMidpointDist
            tableArray{itracts,18}=nan;
            %stDevMidpointDist
            tableArray{itracts,19}=nan;
        else
            %endpointDensity1
            tableArray{itracts,11}=results.WBFG.tractStats{itracts}.endpointDensity1;
            %endpointDensity2
            tableArray{itracts,12}=results.WBFG.tractStats{itracts}.endpointDensity2;
            %avgEndpointDist1
            tableArray{itracts,13}=results.WBFG.tractStats{itracts}.avgEndpointDist1;
            %avgEndpointDist2
            tableArray{itracts,14}=results.WBFG.tractStats{itracts}.avgEndpointDist2;
            %stDevEndpointDist1
            tableArray{itracts,15}=results.WBFG.tractStats{itracts}.stDevEndpointDist1;
            %stDevEndpointDist2
            tableArray{itracts,16}=results.WBFG.tractStats{itracts}.stDevEndpointDist2;
            %midpointDensity
            tableArray{itracts,17}=results.WBFG.tractStats{itracts}.midpointDensity;
            %avgMidpointDist
            tableArray{itracts,18}=results.WBFG.tractStats{itracts}.avgMidpointDist;
            %stDevMidpointDist
            tableArray{itracts,19}=results.WBFG.tractStats{itracts}.stDevMidpointDist;
        end
        %norms.volumeProp
        tableArray{itracts,20}=results.WBFG.tractStats{itracts}.norms.volumeProp;
        %norms.countProp
        tableArray{itracts,21}=results.WBFG.tractStats{itracts}.norms.countProp;
        %norms.wireProp
        tableArray{itracts,22}=results.WBFG.tractStats{itracts}.norms.wireProp;
    end
    tableOut = cell2table(tableArray,...
    'VariableNames',fullFieldNames);
end
writetable(tableOut,'output_FiberStats.csv')
end