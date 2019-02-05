function[figHandle, results]= bsc_feAndSegQualityCheck(feORwbfg, classification, saveDir)
%function[figHandle, results]= bsc_feAndAFQqualityCheck(fe, classification, saveDir)
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
%       of your input FE structure and (if included) AFQ segmentation.  See
%       code for specifics of all fields and their interpretation
%
% % (C) Daniel Bullock, 2017, Indiana University
%% preliminaries
rangeSet=11:250;
% loads AFQ classification file if it is a path
if ~notDefined('classification')
    if ischar(classification)
        load(classification)
        %catches classification if it was saved.  Can probably be made more
        %robust to other saving schemas
    end
end

[wbFG, fe] = bsc_LoadAndParseFiberStructure(feORwbfg);

%a test for the
% if isempty(find(isnan(fe.life.fit.weights),1))>0
%     nanNum=length(find(isnan(fe.life.fit.weights)));
%     fprintf('\n %i NaN values detected in fe.life.fit.weights', nanNum)
% end

%if classification exists do the tract analysis, otherwise just do the
%WBFG
if ~notDefined('classification')
    [results]= wma_quantAllWMNorm(feORwbfg,classification);
else
    [results]= wma_quantWBFG(feORwbfg);
end

%% huge set of conditionals for plotting
% if both fe and class are present
figure
if ~notDefined('fe') & ~notDefined('classification')
    %readout stats
    textBoxHandle=subplot(2,5,[1 6]);
    textBoxPos=get(textBoxHandle,'position');
    annotation('textbox',...
        [textBoxPos(1) textBoxPos(2) textBoxPos(3) textBoxPos(4)],...
        'String',{['Raw Stream Count: ', num2str(results.WBFG.stream_count)],...
        ['Total Stream Volume: ', num2str(results.WBFG.volume)],...
        ['Raw fit curve A: ', num2str(results.WBFG.LogFitA)],...
        ['Raw fit curve B: ', num2str(results.WBFG.LogFitB)],...
        ['Raw Mean stream length: ', num2str(results.WBFG.avg_stream_length)],...
        ['Raw stream length stDev: ', num2str(results.WBFG.stream_length_stdev)],...
        ['LiFE RMSE: ', num2str(results.LiFEstats.RMSE.WB)],...
        ['LiFE All Vox Error: ', num2str(results.LiFEstats.RMSE.WB_norm_total)],...
        ['LiFE Survivors count: ', num2str(results.LiFEstats.posWBFG.stream_count)],...
        ['LiFE Survivors pct: ', num2str((results.LiFEstats.posWBFG.stream_count/results.WBFG.stream_count)*100)],...
        ['Nonzero fit curve A: ', num2str(results.LiFEstats.posWBFG.LogFitA)],...
        ['Nonzero fit curve B: ', num2str(results.LiFEstats.posWBFG.LogFitB)],...
        ['Raw Identified Streams: ', num2str(sum(classification.index>0))],...
        ['Raw Identified Streams Proportion: ', num2str(sum(classification.index>0)/results.WBFG.stream_count)]},...
        'FontSize',14,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    %unused
    %['Avg stream Asym ratio: ', num2str(results.WBFG.avgAsymRat)],...
    %['stream Asym ratio stDev: ', num2str(results.WBFG.avgAsymRat)],...
    % ['Tracts Segmented: ', num2str(sum(length(rightNames)+length(interHemiNames)))],...
    
    posBool=fe.life.fit.weights>0;
    
    % Plot comparing pre and post life streamline proportion by length
    subplot(2,5,2)
    hold on
    plot ((rangeSet),results.WBFG.lengthProps(rangeSet),'b', 'LineWidth',1.25)
    plot ((rangeSet),(results.WBFG.LogFitA*exp(results.WBFG.LogFitB*(rangeSet-10))),'c', 'LineWidth',1.25)
    %plot ((rangeSet),lognpdf(rangeSet,results.WBFG.LogFitA,results.WBFG.LogFitB),'c', 'LineWidth',1.25)
    
    plot ((rangeSet),results.LiFEstats.posWBFG.lengthProps(rangeSet),'r', 'LineWidth',1.25)
    plot ((rangeSet),(results.LiFEstats.posWBFG.LogFitA*exp(results.LiFEstats.posWBFG.LogFitB*(rangeSet-10))),'m', 'LineWidth',1.25)
    
    
    title({'Normalized WBFG & NonZero Weighted',' Stream Count Comparison'})
    legend('WBFG','WBFG Fit','NonZero','NonZero Fit')
    xlabel('Streamline Length (mm)')
    ylabel('Whole Brain proportion')
    
    % Plot illustrating the bias in the validated streamlines as assocaited
    % with length
    subplot(2,5,3)
    hold on
    plot ((rangeSet),(results.WBFG.lengthProps(rangeSet)-results.LiFEstats.posWBFG.lengthProps(rangeSet))*10000,'g', 'LineWidth',1.25)
    plot ((rangeSet),(zeros(1,length(rangeSet)))*10000,'r', 'LineWidth',1.25)
    title({'Survival Bias','Relative to Streamline Length'})
    ylim([-50,50])
    legend('WBFG ratio - Validated ratio','No Bias')
    xlabel('Streamline Length (mm)')
    ylabel('Survival bias (%)')
    
    %plot looking at cumulative fiber length proportions
    cumValid=zeros(1,length(results.WBFG.lengthProps));
    cumWBFG=cumValid;
    for ilengths=5:length(results.WBFG.lengthProps)
        cumValid(ilengths)=results.LiFEstats.posWBFG.lengthProps(ilengths)+sum(cumValid(ilengths-1));
        cumWBFG(ilengths)=results.WBFG.lengthProps(ilengths)+sum(cumWBFG(ilengths-1));
    end
    
    subplot(2,5,7)
    hold on
    plot (((rangeSet)),cumWBFG((rangeSet)),'b', 'LineWidth',1.25)
    plot (((rangeSet)),cumValid((rangeSet)),'r', 'LineWidth',1.25)
    title({'Cumulative portion of fibers','in connectome, by length'})
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Portion of tracts less than or equal to length')
    
    %% Seg Case -- Plots Specific to segmentation output
    
    %set values
    classBool=classification.index>0;
    [WBFGclassHist, ~]=histcounts(results.WBFG.lengthData(classBool),(1:300));
    [validClassHist, ~]=histcounts(results.WBFG.lengthData(classBool & posBool),(1:300));
    
    % Plot comparing pre and post life classified streamline proportion by length
    subplot(2,5,8)
    hold on
    plot ((rangeSet),(WBFGclassHist(rangeSet)/results.WBFG.stream_count)*100,'b', 'LineWidth',1.25)
    plot ((rangeSet),(validClassHist(rangeSet)/results.LiFEstats.posWBFG.stream_count)*100,'r', 'LineWidth',1.25)
    title({'Classified Streamline Proportion','Comparison: WBFG & Surviving'})
    legend('WBFG, AFQ classified','Validated & AFQ classified')
    xlabel('Streamline Length (mm)')
    ylabel('Proportion of Whole Brain Streamlines Classified (%)')
    
    
    
    %% computation for seg bar plots
    %count plot
    classificationGrouping = wma_classificationStrucGrouping(classification);
    countPlotInput=zeros(2,(length(classificationGrouping.names)));
    volPlotInput=zeros(2,(length(classificationGrouping.names)));
    
    for itracts=1:length(classificationGrouping.names)
        fprintf('\n %s',classificationGrouping.names{itracts})
        tractIndexes=unique(classification.index(classificationGrouping.index==itracts));
        for iVariants=1:length(tractIndexes)
            countPlotInput(iVariants,itracts)=results.WBFG.tractStats{tractIndexes(iVariants)}.norms.countProp;
            volPlotInput(iVariants,itracts)=results.WBFG.tractStats{tractIndexes(iVariants)}.norms.volumeProp;
        end
    end
    
    bottomPlot1=subplot(2,5,[4,5]);
    hold on
    bar((countPlotInput')*100)
    title({'Proportion of connectome',' streamlines in tract'})
    legend('Left','Right')
    xlabel('Tract')
    ylabel('% classificaiton input streamlines in tract (%)')
    ylim([-0 max(max(countPlotInput))*133])
    set(gca,'xtick',[1:1:length(classificationGrouping.names)])
    set(gca,'XTickLabel',classificationGrouping.names, 'FontSize',8,'FontName','Times')
    bottomPlot1.XTickLabelRotation=-45;
    
    bottomPlot2=subplot(2,5,[9,10]);
    hold on
    bar((volPlotInput')*100)
    title({'Proportion of wm volume','occupied by tract'})
    legend('Left','Right')
    xlabel('Tract')
    ylabel('% wm volume proportion occupied by tract (%)')
    ylim([-0 max(max(volPlotInput))*133])
    set(gca,'xtick',[1:1:length(classificationGrouping.names)])
    set(gca,'XTickLabel',classificationGrouping.names, 'FontSize',8,'FontName','Times')
    bottomPlot2.XTickLabelRotation=-45;
    
    fig = gcf;
    fig.Position = [100 50 2200 1250];
    
    % if only fe is passed
elseif ~notDefined('fe') & notDefined('classification')
    textBoxHandle=subplot(3,2,[1 3 5]);
    textBoxPos=get(textBoxHandle,'position');
    annotation('textbox',...
        [textBoxPos(1) textBoxPos(2) textBoxPos(3) textBoxPos(4)],...
        'String',{['Raw Stream Count: ', num2str(results.WBFG.stream_count)],...
        ['Total Stream Volume: ', num2str(results.WBFG.volume)],...
        ['Raw fit curve A: ', num2str(results.WBFG.LogFitA)],...
        ['Raw fit curve B: ', num2str(results.WBFG.LogFitB)],...
        ['Raw Mean stream length: ', num2str(results.WBFG.avg_stream_length)],...
        ['Raw stream length stDev: ', num2str(results.WBFG.stream_length_stdev)],...
        ['LiFE RMSE: ', num2str(results.LiFEstats.RMSE.WB)],...
        ['LiFE All Vox Error: ', num2str(results.LiFEstats.RMSE.WB_norm_total)],...
        ['Nonzero fit curve A: ', num2str(results.LiFEstats.posWBFG.LogFitA)],...
        ['Nonzero fit curve B: ', num2str(results.LiFEstats.posWBFG.LogFitB)],...
        ['LiFE Survivors count: ', num2str(results.LiFEstats.posWBFG.stream_count)],...
        ['LiFE Survivors pct: ', num2str((results.LiFEstats.posWBFG.stream_count/results.WBFG.stream_count)*100)]},...
        'FontSize',14,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    %unused
    %['Avg stream Asym ratio: ', num2str(results.WBFG.avgAsymRat)],...
    %['stream Asym ratio stDev: ', num2str(results.WBFG.avgAsymRat)],...
    % ['Tracts Segmented: ', num2str(sum(length(rightNames)+length(interHemiNames)))],...
    
    posBool=fe.life.fit.weights>0;
    
    % Plot comparing pre and post life streamline proportion by length
    subplot(3,2,2)
    hold on
    plot ((rangeSet),results.WBFG.lengthProps(rangeSet),'b', 'LineWidth',1.25)
    plot ((rangeSet),(results.WBFG.LogFitA*exp(results.WBFG.LogFitB*(rangeSet-10))),'c', 'LineWidth',1.25)
    %plot ((rangeSet),lognpdf(rangeSet,results.WBFG.LogFitA,results.WBFG.LogFitB),'c', 'LineWidth',1.25)
    
    plot ((rangeSet),results.LiFEstats.posWBFG.lengthProps(rangeSet),'r', 'LineWidth',1.25)
    plot ((rangeSet),(results.LiFEstats.posWBFG.LogFitA*exp(results.LiFEstats.posWBFG.LogFitB*(rangeSet-10))),'m', 'LineWidth',1.25)
    
    
    title({'Normalized WBFG & NonZero Weighted',' Stream Count Comparison'})
    legend('WBFG','WBFG Fit','NonZero','NonZero Fit')
    xlabel('Streamline Length (mm)')
    ylabel('Whole Brain proportion')
    
    % Plot illustrating the bias in the validated streamlines as assocaited
    % with length
    subplot(3,2,4)
    hold on
    plot ((rangeSet),(results.WBFG.lengthProps(rangeSet)-results.LiFEstats.posWBFG.lengthProps(rangeSet))*10000,'g', 'LineWidth',1.25)
    plot ((rangeSet),(zeros(1,length(rangeSet)))*10000,'r', 'LineWidth',1.25)
    title({'Survival Bias','Relative to Streamline Length'})
    ylim([-50,50])
    legend('WBFG ratio - Validated ratio','No Bias')
    xlabel('Streamline Length (mm)')
    ylabel('Survival bias (%)')
    
    %plot looking at cumulative fiber length proportions
    cumValid=zeros(1,length(results.WBFG.lengthProps));
    cumWBFG=cumValid;
    for ilengths=5:length(results.WBFG.lengthProps)
        cumValid(ilengths)=results.LiFEstats.posWBFG.lengthProps(ilengths)+sum(cumValid(ilengths-1));
        cumWBFG(ilengths)=results.WBFG.lengthProps(ilengths)+sum(cumWBFG(ilengths-1));
    end
    
    subplot(3,2,6)
    hold on
    plot (((rangeSet)),cumWBFG((rangeSet)),'b', 'LineWidth',1.25)
    plot (((rangeSet)),cumValid((rangeSet)),'r', 'LineWidth',1.25)
    title({'Cumulative portion of fibers','in connectome, by length'})
    legend('WBFG','Validated')
    xlabel('Streamline Length (mm)')
    ylabel('Portion of tracts less than or equal to length')
    
    fig = gcf;
    fig.Position = [100 100 1200 1000];
    
    %% no life, but with classification
elseif notDefined('fe') & ~notDefined('classification')
    %readout stats
    textBoxHandle=subplot(3,4,[1 5 9]);
    textBoxPos=get(textBoxHandle,'position');
    annotation('textbox',...
        [textBoxPos(1) textBoxPos(2) textBoxPos(3) textBoxPos(4)],...
        'String',{['Raw Stream Count: ', num2str(results.WBFG.stream_count)],...
        ['Total Stream Volume: ', num2str(results.WBFG.volume)],...
        ['Raw fit curve A: ', num2str(results.WBFG.LogFitA)],...
        ['Raw fit curve B: ', num2str(results.WBFG.LogFitB)],...
        ['Raw Mean stream length: ', num2str(results.WBFG.avg_stream_length)],...
        ['Raw stream length stDev: ', num2str(results.WBFG.stream_length_stdev)],...
        ['Raw Identified Streams: ', num2str(sum(classification.index>0))],...
        ['Raw Identified Streams Proportion: ', num2str(sum(classification.index>0)/results.WBFG.stream_count)]},...
        'FontSize',14,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    %unused
    %['Avg stream Asym ratio: ', num2str(results.WBFG.avgAsymRat)],...
    %['stream Asym ratio stDev: ', num2str(results.WBFG.avgAsymRat)],...
    % ['Tracts Segmented: ', num2str(sum(length(rightNames)+length(interHemiNames)))],...
    
    %posBool=fe.life.fit.weights>0;
    
    % Plot comparing pre and post life streamline proportion by length
    subplot(3,4,2)
    hold on
    plot ((rangeSet),results.WBFG.lengthProps(rangeSet),'b', 'LineWidth',1.25)
    plot ((rangeSet),(results.WBFG.LogFitA*exp(results.WBFG.LogFitB*(rangeSet-10))),'c', 'LineWidth',1.25)
    %plot ((rangeSet),lognpdf(rangeSet,results.WBFG.LogFitA,results.WBFG.LogFitB),'c', 'LineWidth',1.25)
    
    title({'Normalized WBFG & NonZero Weighted',' Stream Count Comparison'})
    legend('WBFG','WBFG Fit')
    xlabel('Streamline Length (mm)')
    ylabel('Whole Brain proportion')
    
    %plot looking at cumulative fiber length proportions
    cumValid=zeros(1,length(results.WBFG.lengthProps));
    cumWBFG=cumValid;
    for ilengths=5:length(results.WBFG.lengthProps)
        cumWBFG(ilengths)=results.WBFG.lengthProps(ilengths)+sum(cumWBFG(ilengths-1));
    end
    
    subplot(3,4,3)
    hold on
    plot (((rangeSet)),cumWBFG((rangeSet)),'b', 'LineWidth',1.25)
    title({'Cumulative portion of fibers','in connectome, by length'})
    legend('WBFG')
    xlabel('Streamline Length (mm)')
    ylabel('Portion of tracts less than or equal to length')
    
    %% Seg Case -- Plots Specific to segmentation output
    
    %set values
    classBool=classification.index>0;
    [WBFGclassHist, ~]=histcounts(results.WBFG.lengthData(classBool),(1:300));
    
    % Plot comparing pre and post life classified streamline proportion by length
    subplot(3,4,4)
    hold on
    plot ((rangeSet),(WBFGclassHist(rangeSet)/results.WBFG.stream_count)*100,'b', 'LineWidth',1.25)
    title({'Classified Streamline Proportion','Comparison: WBFG & Surviving'})
    legend('WBFG, classified')
    xlabel('Streamline Length (mm)')
    ylabel('Proportion of Whole Brain Streamlines Classified (%)')
    
    
    
    %% computation for seg bar plots
    %count plot
    classificationGrouping = wma_classificationStrucGrouping(classification);
    countPlotInput=zeros(2,(length(classificationGrouping.names)));
    volPlotInput=zeros(2,(length(classificationGrouping.names)));
    
    for itracts=1:length(classificationGrouping.names)
        fprintf('\n %s',classificationGrouping.names{itracts})
        tractIndexes=unique(classification.index(classificationGrouping.index==itracts));
        for iVariants=1:length(tractIndexes)
            countPlotInput(iVariants,itracts)=results.WBFG.tractStats{tractIndexes(iVariants)}.norms.countProp;
            volPlotInput(iVariants,itracts)=results.WBFG.tractStats{tractIndexes(iVariants)}.norms.volumeProp;
        end
    end
    
    bottomPlot1=subplot(3,4,[6 7 8]);
    hold on
    bar((countPlotInput')*100)
    title({'Proportion of connectome',' streamlines in tract'})
    legend('Left','Right')
    xlabel('Tract')
    ylabel('% classificaiton input streamlines in tract (%)')
    ylim([-0 max(max(countPlotInput))*133])
    set(gca,'xtick',[1:1:length(classificationGrouping.names)])
    set(gca,'XTickLabel',classificationGrouping.names, 'FontSize',8,'FontName','Times')
    bottomPlot1.XTickLabelRotation=-45;
    
    bottomPlot2=subplot(3,4,[10 11 12]);
    hold on
    bar((volPlotInput')*100)
    title({'Proportion of wm volume','occupied by tract'})
    legend('Left','Right')
    xlabel('Tract')
    ylabel('% wm volume proportion occupied by tract (%)')
    ylim([-0 max(max(volPlotInput))*133])
    set(gca,'xtick',[1:1:length(classificationGrouping.names)])
    set(gca,'XTickLabel',classificationGrouping.names, 'FontSize',8,'FontName','Times')
    bottomPlot2.XTickLabelRotation=-45;
    
    fig = gcf;
    fig.Position = [100 100 1600 1000];
    %% no classification no fe
elseif notDefined('fe') & notDefined('classification')
    %readout stats
    textBoxHandle=subplot(2,2,[1 3]);
    textBoxPos=get(textBoxHandle,'position');
    annotation('textbox',...
        [textBoxPos(1) textBoxPos(2) textBoxPos(3) textBoxPos(4)],...
        'String',{['Raw Stream Count: ', num2str(results.WBFG.stream_count)],...
        ['Total Stream Volume: ', num2str(results.WBFG.volume)],...
        ['Raw fit curve A: ', num2str(results.WBFG.LogFitA)],...
        ['Raw fit curve B: ', num2str(results.WBFG.LogFitB)],...
        ['Raw Mean stream length: ', num2str(results.WBFG.avg_stream_length)],...
        ['Raw stream length stDev: ', num2str(results.WBFG.stream_length_stdev)]},...
        'FontSize',14,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 0],...
        'LineWidth',2,...
        'BackgroundColor',[0.9  0.9 0.9],...
        'Color',[0.84 0.16 0]);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    
    %unused
    %['Avg stream Asym ratio: ', num2str(results.WBFG.avgAsymRat)],...
    %['stream Asym ratio stDev: ', num2str(results.WBFG.avgAsymRat)],...
    % ['Tracts Segmented: ', num2str(sum(length(rightNames)+length(interHemiNames)))],...
    
    %posBool=fe.life.fit.weights>0;
    
    % Plot comparing pre and post life streamline proportion by length
    subplot(2,2,2)
    hold on
    plot ((rangeSet),results.WBFG.lengthProps(rangeSet),'b', 'LineWidth',1.25)
    plot ((rangeSet),(results.WBFG.LogFitA*exp(results.WBFG.LogFitB*(rangeSet-10))),'c', 'LineWidth',1.25)
    %plot ((rangeSet),lognpdf(rangeSet,results.WBFG.LogFitA,results.WBFG.LogFitB),'c', 'LineWidth',1.25)
    
    title({'Normalized WBFG',' Stream Count Comparison'})
    legend('WBFG','WBFG Fit')
    xlabel('Streamline Length (mm)')
    ylabel('Whole Brain proportion')
    
    %plot looking at cumulative fiber length proportions
    cumValid=zeros(1,length(results.WBFG.lengthProps));
    cumWBFG=cumValid;
    for ilengths=5:length(results.WBFG.lengthProps)
        cumWBFG(ilengths)=results.WBFG.lengthProps(ilengths)+sum(cumWBFG(ilengths-1));
    end
    
    subplot(2,2,4)
    hold on
    plot (((rangeSet)),cumWBFG((rangeSet)),'b', 'LineWidth',1.25)
    title({'Cumulative portion of fibers','in connectome, by length'})
    legend('WBFG')
    xlabel('Streamline Length (mm)')
    ylabel('Portion of tracts less than or equal to length')
    
    fig = gcf;
    fig.Position = [100 100 1400 600];
end

%gets figure handle
figHandle=gcf;

if ~notDefined('saveDir')
    saveas(gcf,strcat(saveDir,'/tractomeStatPlot.epsc'));
    save(strcat(saveDir,'/tractomeResultStruc.mat'),'results')
end

end