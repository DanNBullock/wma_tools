function callMatrix = wma_makeGitDirTree(targetDir)

wma_makeGitUtilityDir(targetDir)

sourceList=dir(targetDir);
dirBool=[sourceList.isdir];
nameList={sourceList(:).name};
for inames=1:length(nameList)
    dotBool(inames)= strcmp(nameList{inames}(1),'.');
end

dirIndexes=find(dirBool&~dotBool);

allFiles=[];

for iFiles=1:length(dirIndexes)
    dirPath=fullfile(targetDir,sourceList(dirIndexes(iFiles)).name);
    dirContents=dir(dirPath);

    ignoreBool=or(cell2mat(cellfun(@(x) strcmp(x,'.'),{dirContents(:).name},'UniformOutput',false)),...
    cell2mat(cellfun(@(x) strcmp(x,'..'),{dirContents(:).name},'UniformOutput',false)));
    allFiles=horzcat(allFiles, {dirContents(~ignoreBool).name});
end


callMatrix=zeros(length(allFiles));
%% Matrix construction
%row=source, column = target
for iUniqueFiles=1:length(allFiles)
    
    dependancyList=matlab.codetools.requiredFilesAndProducts(allFiles{iUniqueFiles},'toponly'); 
   
    dependancyNames=[];
    
    for iDependancies=1:length(dependancyList)
        [aa bb cc]=fileparts(dependancyList{iDependancies});
        dependancyNames=horzcat(dependancyNames,{strcat(bb,cc)});
    end
    
    boolVec=[];
    
    for iFiles=1:length(allFiles)
        boolVec(iFiles)=sum(cell2mat(cellfun(@(x) strcmp(x,allFiles{iFiles}),dependancyNames,'UniformOutput',false)));
    end
    
    callMatrix(iUniqueFiles,:)=boolVec;
   
    fprintf('\n %i Functon dependancies of %s determined',length(dependancyList),allFiles{iUniqueFiles})
end

figure
dg=digraph(callMatrix, allFiles,'OmitSelfLoops');
plot(dg)
inputFunctions = indegree(dg);
figure
etreeplot(callMatrix)


end

