function wma_makeFunctionListTxt(functionName)


if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
end

functionPath=which(functionName);


[filepath,name,ext] = fileparts(functionPath);
sourceDir=filepath;
dirName=name;

dependancyList=matlab.codetools.requiredFilesAndProducts(functionPath);

sourceList=dir(sourceDir);
fileBool=~[sourceList.isdir];
fileIndexes=find(fileBool);

for iFiles=1:length(fileIndexes)
    fileName{iFiles}=fullfile(sourceDir,sourceList(fileIndexes(iFiles)).name);
end

intersection = fileName(ismember(fileName, dependancyList));

callMatrix=zeros(length(intersection));
%% Matrix construction
%row=source, column = target
for iUniqueFiles=1:length(intersection)
    
    [FILEPATH,NAME,EXT]=fileparts(intersection{iUniqueFiles});
    NamesVec{iUniqueFiles}=NAME;
    
    dependancyList=matlab.codetools.requiredFilesAndProducts(intersection{iUniqueFiles},'toponly'); 
   
    dependancyNames=[];
    
    
    boolVec=[];
    
    for iFiles=1:length(intersection)
        boolVec(iFiles)=sum(cell2mat(cellfun(@(x) strcmp(x,intersection{iFiles}),dependancyList,'UniformOutput',false)));
    end
    
    callMatrix(iUniqueFiles,:)=boolVec;
   
    fprintf('\n %i Functon dependancies of %s determined',length(dependancyList),intersection{iUniqueFiles})
end

%figure('visible', 'off')
dg=digraph(callMatrix, NamesVec,'OmitSelfLoops');
plot(dg)
% inputFunctions = indegree(dg);
% figure
% etreeplot(callMatrix)

saveas(gcf,strcat(pwd,strcat(dirName,'map.png')));

T = cell2table(intersection');
 
% Write the table to a CSV file
writetable(T,strcat(pwd,strcat(dirName,'functionCalls.csv')));

%close all

end



