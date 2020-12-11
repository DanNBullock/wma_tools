function [dataDirPaths]=bsc_getBrainlifeDataDirPaths(projectDir,searchStrings)
%dataDirPaths]=bsc_getBrainlifeDataDirPaths(projectDir,searchStrings)
%
% This function looks across a project directory and obtains the sub-paths
% to the relevant data directories across subjects.
%
%
%  Inptus
%
% projectDir:  path to the top level data directory for the project
%
% searchStrings:  a cell array of strings which must ALL be present in the
% directory name for the data modality you wish to generate paths to.  Can
% include multiple strings, for example the data modality and any relevant
% tags
%
% Output
%
% dataDirPaths:  A cell array of paths to the relevant data directories
%
%  Written By Daniel Bullock June 7 2020
%% Begin code

%get the project dir contents, should be a list of subjects and other stuff
projectDirContents=dir(projectDir);

%get a list of the names
contentNames={projectDirContents.name};

%extract only the subject data directories, in case any other directories
%are in there
subjectDirNames=contentNames(and(contains(contentNames,'sub-'),[projectDirContents.isdir]));

%create output vector
dataDirPaths=cell(length(subjectDirNames),1);

%iterate across subject directories
for iSubjectDirs=1:length(dataDirPaths)
    %build path to current subject directory
    currentSubjDir=fullfile(projectDir,subjectDirNames{iSubjectDirs});
    %directory contents for this subject
    currentSubjDirContent=dir(currentSubjDir);
    %get the names
    currentSubjDirFileNames={currentSubjDirContent.name};
    %create a blank boolean vector for all the requirements
    currentCriteria=false(length(currentSubjDirFileNames),length(searchStrings),1);
    %iterate across the criteria
    for iCriteria=1:length(searchStrings)
        currentCriteria(:,iCriteria)=contains(currentSubjDirFileNames,searchStrings{iCriteria})';
    end
    %now append the isDir requirement
    currentCriteria=horzcat(currentCriteria,[currentSubjDirContent.isdir]');
    %now find the index of the relevant file
    curentDirectoryIndex=find(all(currentCriteria,2));
    %throw an error if multiple are detected
    if length(curentDirectoryIndex)>1
        error('multiple potential data directories found in %s',currentSubjDir)
    elseif isempty(curentDirectoryIndex)
        warning('no valid data directory found in %s',currentSubjDir);
        dataDirPaths{iSubjectDirs}=[];
    else    
    %otherwise use this to form the directory path
    dataDirPaths{iSubjectDirs}=fullfile(currentSubjDir,currentSubjDirFileNames{curentDirectoryIndex});
    end
end

%get the total number of paths created
pathCreatedCount=sum(cellfun(@ischar,dataDirPaths));

fprintf('\n %i valid data paths created for project %s',pathCreatedCount,projectDir)
end  