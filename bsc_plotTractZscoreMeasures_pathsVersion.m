function  bsc_plotTractZscoreMeasures_pathsVersion(csvPaths,plotProperties,saveDir)


%workingDir='/N/dc2/projects/lifebid/HCP/Dan/EcogProject/proj-5c33a141836af601cc85858d'
%identifierTag='cleaned'
%plotProperties=[19 20]

mkdir(fullfile(saveDir,'data'));

%csvPaths = tractStatNamesGen(workingDir,identifierTag)

[avgTable, stdTable]=bsc_tableAverages(csvPaths);

[domainNames,propertyNames,valueArray]= bsc_normalizeStatMeasures_pathsVersion(csvPaths);

%nonNameProperties=propertyNames(2:end);

%onlysubjNames=erase(subjNames,'sub-');

mkdir(fullfile(workingDir,'image/'));

for iplotProperties=1:length(plotProperties)
    leftLabels=[];
    for iDomains=1:length(domainNames)
    leftLabels{iDomains}=[domainNames{iDomains},' (',num2str(avgTable{iDomains,plotProperties(iplotProperties)+1}),' +/- ',num2str(stdTable{iDomains,plotProperties(iplotProperties)+1})];
    end
    figure
    %indexing at 2 to squeeze out wbfg, casuses problems otherwise
    plotArray=squeeze(valueArray(:,plotProperties(iplotProperties)+1,:));
    %plotArray(isnan(plotArray))=0;
    minMax=[abs(min(min(min(plotArray,[],'omitnan')))),max(max(max(plotArray,[],'omitnan')))]
    maxOrMin=max(minMax);
    imagesc(plotArray)
    redBlueMap=redblue(99999);
    %hyper inelegant, relies on precision to give you black. In some cases,
    %a value sufficiently close to 0 may appear to be zero, and thus nan on
    %this scale.
    %redBlueMap(50000,:)=[0,0,0];

    colormap(redBlueMap);
    caxis([-maxOrMin maxOrMin])
    colormap(redBlueMap);
    set(gca,'TickLength',[0,0])
    yticks([0:1:length(domainNames)])
    %WHY DO I HAVE TO CIRCULAR SHIFT?
    set(gca,'YTickLabel',circshift(leftLabels,1))
    set(gca,'XTickLabel',[]);
    ax = gca;
    ax.YLabel.String = 'Tracts';
     ax.XLabel.String = 'Subjects';
    ax.Title.String = [identifierTag,' tracts'' ', propertyNames{plotProperties(iplotProperties)+1},'s'];
    curMap=colorbar;
    curMap.Label.String = 'Z score';
    figName=[identifierTag,'_', propertyNames{plotProperties(iplotProperties)+1}];
    
    saveas(gcf,[fullfile(workingDir,'image/'),figName,'.svg']);
end

if saveFlag
    save([fullfile(saveDir,'data'),'groupZscoreData.mat'],'subjNames','domainNames','propertyNames','valueArray')
end

end