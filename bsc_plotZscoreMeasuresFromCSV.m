function  bsc_plotZscoreMeasuresFromCSV(csvPaths,plotProperties,saveDir,subSelect)
% function  bsc_plotTractZscoreMeasures_pathsVersion_v2(csvPaths,plotProperties,saveDir)
%
%  Computes zscore for specified properties and plots them in a matrix,
%  with rows corresponding to domains and columns corresponding to
%  subjects (or whatever csvPaths means).
%
%  INPUTS
%  csvPaths:  paths to the csv files which are to be amalgamated into a
%  single structure.
%
%  plotProperties:  either a string corresponding to the exact property to
%  be plotted, or a number indicating its placement IN THE OUTDATASTRUC
%  (this means the number should ignore the first column if that is used
%  for domain labeling)
%
%  saveDir:  the directory to save down all relevant inputs
%
%  subSelect:  a numeric vector indicating which domains are to be selected
%  for plotting.  If not defined, will plot all.
%
%  OUTPUTS
%  NONE - SAVES RELEVANT OUTPUT TO SPECIFIED DIRECTORY
%%  begin code
% create 
mkdir(fullfile(saveDir,'data/'));
mkdir(fullfile(saveDir,'image/'));

%generate compiled 3d data structure from tables
[allProperties, allDomains, outDataStruc]=bsc_csvTables2DataStruc(csvPaths);

%NOTE THIS ASSUMES NORMAL DISTRIBUTION
%compute mean of structure
avgTable = mean(outDataStruc,3,'omitnan');
%compute std of structure
stdTable = std(outDataStruc,0,3,'omitnan');

% get size of output structure
dataStrucDim=size(outDataStruc);
%initialize structure for normalized values
normalStruc=zeros(dataStrucDim);
%iterate through subjects and apply z score transform
for isubjects=1:dataStrucDim(3)
    normalStruc(:,:,isubjects)=rdivide(minus(outDataStruc(:,:,isubjects),avgTable),stdTable);  
end


%convert plot properties to something that the plotting function can use
for iplotProperties=1:length(plotProperties)
    if ischar(plotProperties{iplotProperties})
    plotProperties{iplotProperties}=find(strcmp(plotProperties{iplotProperties},allProperties));
    else
        %probably will error for numbers
        plotProperties{iplotProperties}=plotProperties{iplotProperties};
    end
end

for iplotProperties=1:length(plotProperties)
    leftLabels=[];
    for iDomains=1:length(allDomains)
    leftLabels{iDomains}=[allDomains{iDomains},' (',num2str(avgTable(iDomains,plotProperties{iplotProperties})),' +/- ',num2str(stdTable(iDomains,plotProperties{iplotProperties})),')'];
    end
    figure
    prePlotArray=squeeze(normalStruc(:,plotProperties{iplotProperties},:));
    
    %apply subselect
    if exist('subSelect')
        plotArray=prePlotArray(subSelect,:);
    else
        plotArray=prePlotArray(1:end,:);
    end

    minMax=[abs(min(min(min(plotArray,[],'omitnan')))),max(max(max(plotArray,[],'omitnan')))];
    maxOrMin=max(minMax);
    imagesc(plotArray)
    redBlueMap=redblue(99999);

    colormap(redBlueMap);
    caxis([-maxOrMin maxOrMin])
    colormap(redBlueMap);
    set(gca,'TickLength',[0,0])
    yticks([0:1:length(allDomains)])
    set(gcf, 'Position',  [0, 0, 400+dataStrucDim(3)*3, 300+dataStrucDim(1)*10])
    %WHY DO I HAVE TO CIRCULAR SHIFT?  Try this without for now
    %pause(1)
    set(gca,'YTickLabel',circshift(leftLabels,1))
    %pause(1)
    %set(gca,'YTickLabel',leftLabels)
    set(gca,'XTickLabel',[]);
    ax = gca;
    ax.YLabel.String = 'Tracts';
    ax.XLabel.String = 'Subjects';
    ax.Title.String = [allProperties{plotProperties{iplotProperties}}];
    curMap=colorbar;
    curMap.Label.String = 'Z score';
    figName=[allProperties{plotProperties{iplotProperties}}];
   
    saveas(gcf,[fullfile(saveDir,'image/'),figName,'.svg']);
end

if ~notDefined(saveDir)
    save([fullfile(saveDir,'data/'),'groupZscoreData.mat'],'normalStruc');
    save([fullfile(saveDir,'data'),'rawDataStruc.mat'],'outDataStruc');
end

end