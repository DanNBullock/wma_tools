function[results]= wma_quantWBFG(feORwbfg)
%function[results]= wma_quantWBFG(fe, classification, saveDir)
%
% This function computes a number of statistics for an input fe structure,
% and, if included, an afq segmentation.  if you include a saveDir variable
% it saves down a figure illustrating a number of these statistics
% (relative to length) as well as the results structure.  "Validated" or
% "post-life" refers to streamlines which have been validated by life and
% thus have nonzero weights in fe.life.fit.weights.  Correspondingly,
% "WBFG" or "pre-life" referrs to ALL streamlines from the input whole
% brain tractography.
%
% INPUTS: -feORwbfg:  either a whole brain fiber group / tractogram
% path/object or an fe path / object.
%
% -classification: Either the path to structure or the structure itself.
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
% -saveDir: path that you would like the results struture and plot saved
%       to.  If not defined, does not save the output.
%
% OUTPUTS:
%
% -figHandle: handle for the mutiplot figure.  Can be used to adjust figure
%       properties or to save figure.
%
% -results: a structure with multiple fields summarizing various properties
%       of your input FE structure and (if included) Seg segmentation.  See
%       code for specifics of all fields and their interpretation
%
% % (C) Daniel Bullock, 2017, Indiana University
%% preliminaries

% loads Seg classification file if it is a path
[wbFG, fe] = bsc_LoadAndParseFiberStructure(feORwbfg);

if ~isempty(fe)
    feFlag=true;
else
    feFlag=false;
end

%% resultStrucStuff
% computes length  of positively weighted streamlines
[WBtractStat]= wma_quantTract(wbFG);
results.WBFG=WBtractStat;
results.WBFG.lengthProps=results.WBFG.lengthCounts/results.WBFG.stream_count;

%this may make the structure too large, if so consider removing
streamLengths=zeros(1,length(wbFG.fibers));
for istreamlines=1:length(wbFG.fibers)
    streamLengths(istreamlines)=sum(sqrt(sum(diff(wbFG.fibers{istreamlines},1,2).^2)));
end

results.WBFG.lengthData=streamLengths;

%get stats for lognormal fit
%[paramsOutWB]=lognfit(results.WBFG.lengthProps(results.WBFG.lengthCounts>0));
paramsOutWB = fit((1:length(results.WBFG.lengthProps(results.WBFG.lengthCounts>0)))',results.WBFG.lengthProps(results.WBFG.lengthCounts>0)','exp1');
results.WBFG.LogFitA=paramsOutWB.a;
results.WBFG.LogFitB=paramsOutWB.b;

%% LiFE dependent Analyses
switch feFlag
    case true
        % gets positively weighted streamlines and their indexes
        posIndexes=find(fe.life.fit.weights>0);
        
        %a test for the
        if isempty(find(isnan(fe.life.fit.weights),1))>0
            nanNum=length(find(isnan(fe.life.fit.weights)));
            fprintf('\n %i NaN values detected in fe.life.fit.weights', nanNum)
        end
        
        %get name
        results.LiFEstats.fe_name=fe.name;
        
        %get name
        posWB=wbFG;
        posWB.fibers=wbFG.fibers(posIndexes);
        
        %get tract quantification for
        [posWBtractStat]= wma_quantTract(posWB);
        
        % Computes the average and standard deviation of the voxelwise error. See
        % feGet: 'voxrmses0norm' for more information.
        % Concern: does this include non tractography occupied voxels?
        allvoxelsErr= feGet(fe, 'voxrmses0norm');
        results.LiFEstats.RMSE.voxel_average_=mean(allvoxelsErr(~isnan(allvoxelsErr)));
        results.LiFEstats.RMSE.voxel_stdev=std(allvoxelsErr(~isnan(allvoxelsErr)));
        results.LiFEstats.RMSE.voxel_count=length(allvoxelsErr);
        % Gets the total RMSE for the tractography model.  Not equivalent to
        % sum(allvoxelsErr(~isnan(allvoxelsErr))), hence the above concern.
        results.LiFEstats.RMSE.WB=feGet(fe,'rmsetotal');
        % Whole model error sum
        results.LiFEstats.RMSE.WB_norm_total=sum(allvoxelsErr(~isnan(allvoxelsErr)));
        
        % Computes total length of pre and post LiFE streamline connectome.
        % Compares the values to get proportional reduction in total connectome
        % streamline length.
        results.LiFEstats.posWBFG=posWBtractStat;
        
        results.LiFEstats.compare.survive=results.LiFEstats.posWBFG.stream_count/results.WBFG.stream_count;
        results.LiFEstats.compare.totalWireSurvive=results.LiFEstats.posWBFG.length_total/results.WBFG.length_total;
        %might mess up in cases of zero
        results.LiFEstats.posWBFG.lengthProps=rdivide(results.LiFEstats.posWBFG.lengthCounts,results.LiFEstats.posWBFG.stream_count);
        
        %get stats for lognormal fit
        %paramsOutLiFE=lognfit(posWBtractStat.lengthProps(results.LiFEstats.posWBFG.lengthCounts>0));
        paramsOutLiFE = fit((1:length(results.LiFEstats.posWBFG.lengthProps(results.LiFEstats.posWBFG.lengthCounts>0)))',results.LiFEstats.posWBFG.lengthProps(results.LiFEstats.posWBFG.lengthCounts>0)','exp1');

        results.LiFEstats.posWBFG.LogFitA=paramsOutLiFE.a;
        results.LiFEstats.posWBFG.LogFitB=paramsOutLiFE.b;
        
end

end