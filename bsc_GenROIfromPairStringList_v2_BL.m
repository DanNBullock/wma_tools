function bsc_GenROIfromPairStringList_v2_BL()
% bsc_genNiftiROIfromStringList(feORwbfg,atlas,ROIstring, smoothKernel)
%
% Given a string list of rois (in the specified format) loops over
% the list and generates nifti ROI files for each 
%
%  INPUTS:

%  -atlas: path to an atlas or an atlas.
%
%  -ROIstring:  a string list of atlas based roi specifications from the
%  atlas that you would like merged into a single ROI.  eg: '2 56 30 54; 34
%  654 \n 25 45 56; 23 \n 456 34; 35 75'.  Must correspond to the values in
%  the atlas nifti itself.
%
%  OUTPUTS: 
% -classification: 
% 
% -None, saves output down.
%
%  NOTE:  Makes a call to mri_convert, so requires that FreeSurfer be
%  installed and set up properly.
%
%  (C) Daniel Bullock 2018 Bloomington, Indiana
%% Begin code

 if ~isdeployed
     disp('adding paths');
     addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
     addpath(genpath('/N/u/brlife/git/jsonlab'))
     addpath(genpath('/N/u/brlife/git/vistasoft'))
     addpath(genpath('/N/u/brlife/git/wma_tools'))
 end

%config = loadjson('/N/dc2/projects/lifebid/HCP/Dan/GitStoreDir/ROIs2ROIsSegment/config.json');
config = loadjson('config.json')

ROIstring=config.roiPairs;

if isfield(config,'atlas')
smoothKernel=config.smoothKernel;
end

if isfield(config,'atlas')
       atlas=config.atlas
end

if isfield(config,'rois')
      ROIdirIN=config.rois
end
%alternate, becaue I'm a bad programmer/lazy
if isfield(config,'ROI')
      ROIdirIN=config.ROI
end


%% gen ROI
fprintf('Generating ROIs for the following indicies: \n %s',ROIstring);
%just in case
ROIstring=strrep(ROIstring,'\n',newline);
ROIstring=strrep(ROIstring,';',newline);
stringCells = splitlines(ROIstring);
%% ROI DIR SET
%is this right?  If there was a way to flag a line such that it could be
%revisted later, this would be it, as this is where the troublesome ROI
%directory is made.
roiDirPathOut=fullfile(pwd,'rois/');
mkdir(roiDirPathOut)
%the stem that preceeds the roi identifier (identifier = either integer or
%underscore delimited, concatonated string of integers corresponding to the
%7/20/2019 note:  previously there was some discrepancy with whether there
%was a _ which immediately followed the roi stem (such that ROI_xxx vs
%ROIxxx.  This led to incompatabilities between my and brad's
%functions/appps.
roiStem='ROI';

%Just guessing what the right path for this for now
textDir=pwd;
fileID = fopen(fullfile(textDir,'ROIfiles.txt'),'w');

%% ROI creation loop
for iROIs=1:length(stringCells)
    %% robust input name parsing
    %this will fail if there are underscores
    %here we need begin case parsing
    %basically, we can ignore everything before the first number and just
    %make it conditional to there being underscores in the resultant
    %string.  It is safe to assume that noone using a parcellation will
    %type in something with underscores.  This will work with everything
    %internal to brainlife, but will not work with arbitrarily named
    %uploads from users.
    currChars=num2cell(stringCells{iROIs});
    numCells=cellfun(@str2num,currChars,'UniformOutput',false);
    numBool=~cell2mat(cellfun(@isempty,numCells,'UniformOutput',false));
    %ignoring the stem, whatever it is
    startChar=min(find(numBool));
    curROIStringHold=stringCells{iROIs};
    %what if a literal name inludes numbers...
    curROIStringCut=curROIStringHold(startChar:end);
    if all(isnumeric(str2num(curROIStringCut))) && startChar==1
        %subconditional
        
        if isfield(config,'atlas')
            ROInums=str2num(curROIStringCut);
            dirFlag=false;
            atlasFlag=true;
        else
            ROInames=splitlines(strrep(curROIStringHold,' ',newline));
            ROInames= ROInames(~cellfun('isempty',ROInames));
            dirFlag=true;
            atlasFlag=false;
        end
        %Assume literal input, either concat or literal.  Spaces in name
        %obviously dont work as they would indicate separate file names,
        %according to our conventions
    else
        ROInames=splitlines(strrep(curROIStringHold,' ',newline));
        ROInames= ROInames(~cellfun('isempty',ROInames));
        dirFlag=true;
        atlasFlag=false;
    end
    %% run the merge roi function
    %I think exists would work just as well as notDefined here, and could
    %remove another vistasoft dependancy.  Also, on the off chance that an
    %atlas *and* an roidir is passed, this can be made more robust
    if ~notDefined('atlas')&&atlasFlag
        mergedROI =bsc_roiFromAtlasNums(atlas,ROInums, smoothKernel);
        %operating under presumption that roi.name is a thing...
        currROIName=fullfile(pwd,strcat('/',roiDirPathOut,'/',roiStem,mergedROI.name,'.nii.gz'));
        %write file name to text file
        nameToWrite=erase(currROIName,'.nii.gz');
        fprintf(fileID, strcat(nameToWrite,'\n'))
        [~, ~]=dtiRoiNiftiFromMat (mergedROI,atlas,currROIName,1);
    elseif ~notDefined('ROIdirIN')&&dirFlag
        roiDirContents=dir(ROIdirIN);
        roiDirContentsNames={roiDirContents.name};
        %preemptive charCount for dircontentNames
        roiDirContentsLengths=cellfun(@length,roiDirContentsNames);
        %use concat of .nii.gz to elminiate substring issues.
        assumedROIName=strcat(ROInames,'.nii.gz');
        ROIOutName=[];
        for iROInums=1:length(ROInames)
            presumedFileNameIndex=find(contains(roiDirContentsNames,assumedROIName));
            %check for multiple hits, use heuristic to select
            if length(presumedFileNameIndex)>1
                presumedFileNameIndex=presumedFileNameIndex(min(roiDirContentsLengths(presumedFileNameIndex))==roiDirContentsLengths(presumedFileNameIndex));
            elseif isempty(presumedFileNameIndex)
                warning('no match found for requested roi %s',assumedROIName)
            end
            %if it is still greater than 1, just guess.  Better method
            %might be to pick minimal disagreement
            if length(presumedFileNameIndex)>1
                presumedFileNameIndex=  presumedFileNameIndex(1);
            end
            presumedFileName=roiDirContentsNames{presumedFileNameIndex};
            niiPaths{iROInums}=fullfile(ROIdirIN,presumedFileName);
            periodIndex=strfind(roiDirContentsNames{presumedFileNameIndex},'.');
            if iROInums==1
                ROIOutName=strcat(ROIOutName,presumedFileName(1:periodIndex(1)-1));
            else
                ROIOutName=strcat(ROIOutName,'_',presumedFileName(1:periodIndex(1)-1));
            end
            %now that the paths have been generated, merge
            %removing this conditional, roi should be generated and passed out
            %regardless of singleton nature.
            %if length(ROInums)>1
            % this may result in odd names if merged rois are merged again
            mergedROI = niftiMerge(niiPaths, strcat(roiDirPathOut,ROIOutName,'.nii.gz'));        %end
            currROIName=mergedROI.fname;
            %write file name to text file
            nameToWrite=erase(currROIName,'.nii.gz');
            fprintf(fileID, strcat(nameToWrite,'\n'))
            fprintf('\n saving %s',currROIName)
            niftiWrite(mergedROI,currROIName)
            clear niiPaths
        end
    else
        error('Apparent missing file inputs (atlas or roi dir) relative to current roi specification (%s) ',curROIStringHold)
    end
end
    fclose(fileID)
end