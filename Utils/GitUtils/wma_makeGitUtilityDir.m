function wma_makeGitUtilityDir(targetDir)


% combine with treeplot and etree for great fun
utilDir=fullfile(targetDir,'Utilities');



sourceList=dir(targetDir);
dirBool=[sourceList.isdir];
nameList={sourceList(:).name};
for inames=1:length(nameList)
    dotBool(inames)= strcmp(nameList{inames}(1),'.');
end

dirIndexes=find(dirBool&~dotBool);

allFiles=nameList;

for iFiles=1:length(dirIndexes)
    dirPath=fullfile(targetDir,sourceList(dirIndexes(iFiles)).name);
    dirContents=dir(dirPath);
    dirContentStruc(iFiles).names={dirContents(:).name};
    allFiles=horzcat(allFiles, {dirContents(:).name});
end

uniqueFiles=unique(allFiles);
ignoreBool=or(cell2mat(cellfun(@(x) strcmp(x,'.'),uniqueFiles,'UniformOutput',false)),...
    cell2mat(cellfun(@(x) strcmp(x,'..'),uniqueFiles,'UniformOutput',false)));
uniqueFiles={uniqueFiles{~ignoreBool}};

countVec=zeros(1,length(uniqueFiles));
for iFiles=1:length(dirIndexes)
    for iDirFiles=1:length(dirContentStruc(iFiles).names)
        fileIndex=cellfun(@(x) strcmp(x,dirContentStruc(iFiles).names(iDirFiles)),uniqueFiles,'UniformOutput',false);
        countVec(cell2mat(fileIndex))=countVec(cell2mat(fileIndex))+1;
    end
end

duplicationIndexes=find(countVec>1);

skipcount=0;

if ~isempty (duplicationIndexes)
    mkdir(utilDir);
    for iDuplications=1:length(duplicationIndexes)
        if ~strcmp(uniqueFiles{duplicationIndexes(iDuplications)},'.')
            dupFuncPath=which(uniqueFiles{duplicationIndexes(iDuplications)});
            if ~isempty(dupFuncPath)
            copyfile( dupFuncPath, utilDir);
            else
                fprintf('\n file %s skipped',uniqueFiles{duplicationIndexes(iDuplications)})
            end
        else
            skipcount=skipcount+1;
        end
    end
    fprintf('Utility directory created with %i files',length(duplicationIndexes)-skipcount)
else
    warning('No duplicated files found. No Utility folder created');
end
% FIX THIS, IT TIRES TO DELETE OTHER FOUND FILES
for iFiles=1:length(dirIndexes)
    for iDuplications=1:length(duplicationIndexes)
        dirPath=fullfile(targetDir,sourceList(dirIndexes(iFiles)).name);
        if ~strcmp(dirPath,utilDir)
            dirContents=dir(dirPath);
            dirFileNames={dirContents(:).name};
            
            fileFlag=sum(cell2mat(cellfun(@(x) strcmp(x,uniqueFiles{duplicationIndexes(iDuplications)}),dirFileNames,'UniformOutput',false)));
            if fileFlag
            deleteFileName=fullfile(dirPath,uniqueFiles{duplicationIndexes(iDuplications)});   
            delete( deleteFileName)
            fprintf('\n %s deleted',deleteFileName)
            end
        end
    end
end

end
