function bsc_amalgamateROISacrossDirectories(roiDirs,subselect,outdir,instruction,thresh)
%bsc_amalgamateROISacrossDirectories(roiDirs,subselect,outdir,instruction,thresh)
%
%  Purpose:  This function iterates across input roi directories and
%  creates an amalgum across all rois of the same name across input
%  directories.  IMPORTANT NOTE: given that this is a simple amalgamation
%  of rois, it is presumed that all rois of the same name are in the same
%  reference space and are the same size.
%
%  Inputs
%
%  roiDirs: cell structure with the roi directory paths
%
%  subselect:  a cell structure with strings corresponding to the subset of
%  rois that you would like to amalgamate.  This function will iterate
%  across this variable, and will amalgamate all same-named rois that
%  contain this string.  As an example, inputting 'Arc', would amalgamate
%  rois with the string 'pArc' or 'Arc' in them, including left and right
%  variants
%
%  outdir:  The directory you would like to output the amalgamated rois
%  into
%
%  instruction:  how you would like the amalgum made, either 'sum' or
%  'binarized'.  In the case of 'sum' will simply add the contents of the
%  rois together.  In the case of 'binarized' will mask (i.e. set all nonzero
%  entries to 1), and then sum across directories
%
%  thresh:  the threshold applied before either the summing or binarization
%  takes place.  In this fashion, all entries below this value are set to 0
%  before the amalgamation operation takes place.
%
%  Outputs
%  none, saves down output
%
%  Written by Daniel Bullock June 7 2020
%%  begin code

%set the default thresh value
if isempty(thresh)
    thresh=0;
end

%make the out dir if it doesn't exist
if ~isdir(outdir)
    mkdir(outdir)
end

%create a blank names vec for the rois
roiNamesVec=[];

%do a quick loop over input directories to get the unique roi names
%we do this in case there are different sets of rois across the directories
for iInputDirs=1:length(roiDirs)
    currentDirContents=dir(roiDirs{iInputDirs});
    %get the names
    currentNames={currentDirContents.name};
    %subset it down
    currentNames=currentNames(contains(currentNames,'.nii.gz'));
    roiNamesVec=horzcat(roiNamesVec,currentNames);
end

%find the unique names
uniqueRoiNames=unique(roiNamesVec);

%create a vector for the rois to amalgamate
roisToAmalgamate=[];

%if there are subset inputs, find out which rois meet those criteria
if ~isempty(subselect)
    %loop across criteria
    for iCriteria=1:length(subselect)
        %find which names meet this critiera
        currentValidNames=uniqueRoiNames(contains(uniqueRoiNames,subselect{iCriteria}));
        %add it to the output vector
        roisToAmalgamate=horzcat(roisToAmalgamate,currentValidNames);
    end
    %subset it down to the unique entries
    roisToAmalgamate=unique(roisToAmalgamate);
else
    %otherwise just do all the unique ones
    roisToAmalgamate=uniqueRoiNames;
end
    
%loop across rois
for iAmalgums=1:length(roisToAmalgamate)
    %loop across input directories
    currentAmalgum=[];
    for iInputDirs=1:length(roiDirs)
        %make putative path to first ROI
        currentROIpath=fullfile(roiDirs{iInputDirs},roisToAmalgamate{iAmalgums});
        if isfile(currentROIpath)
            currentROI=niftiRead(currentROIpath);
            %apply the threshold
            currentROI.data(currentROI.data<thresh)=0;
            %depending on the input instruction
            switch instruction
                %either take the current roiData as is
                case 'sum'
                    preparedROI=currentROI;
                    
                    %or binarize it
                case 'binarized'
                    currentROI.data(currentROI.data>0)=1;
                    preparedROI=currentROI;
            end
            %if this is the first instance of this roi, and thus there is
            %no current amalgum object, set this prepared roi to the
            %amalgum
            if isempty(currentAmalgum)
                currentAmalgum=preparedROI;
            else
                %otherwise, if an amalgum already exists, sum them
                currentAmalgum.data=currentAmalgum.data+preparedROI.data;
            end
            %the case in which the roi exists is complete
        else %if the roi doesnt exist
            %do nothing, there's no valid roi for this one
        end
    end
    fprintf('\n file %s complete',roisToAmalgamate{iAmalgums})
    currentSavePath=fullfile(outdir,roisToAmalgamate{iAmalgums});
    niftiWrite(currentAmalgum,currentSavePath);
end   
