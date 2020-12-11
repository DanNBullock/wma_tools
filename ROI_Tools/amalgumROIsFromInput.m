function outputROIS=amalgumROIsFromInput(atlasORroiDir,multipleRequestedAmalgums,smoothKernel)
% outputROIS=amalgumROIsFromInput(atlasORroiDir,multipleRequestedAmalgums,smoothKernel)
%
%  This function either extracts amalgum (in that there are multiple labels),
% or singular (in that there is a single label) ros fom an atlas; OR if a
% directory, will pass singular (in that a single roi is requested; i.e.
% simple passthrough) or amalgum (in that multiple rois are merged into
% one) rois out.
%
%  Inputs
%
%  atlasORroiDir:  Path directly to atlas, or path to roi directory, or an
%  atlas itself
%
%  multipleRequestedAmalgums:  A cell array, each cell of which contains:
%  
%  >  If an roi directory is passed in:  a string indicating the requested roi
%  amalgums.  This function will be looking for a sequence of filenames for
%  rois separated by spaces (eg 'filename1 filename2 filename3').  Can have
%  .nii.gz added to end of name or not, function will add if not present.
%
%  > If an atlas is passed in: either a integer numerical vector or a
%  string representing a sequence of integer representations (which will be
%  converted into a numerical vector using str2num).  Either input method
%  should work fine.
%
% OUTPUTS
%
% A cell structure which contains the requested singular or amalgum rois in
% each cell (sequenced in accordance with the input requests
%
% Dan Bullock 4/23/20
%% begin code

%load the atlas if an atlas is passed

if isstruct(atlasORroiDir)
    %maybe implement more checks here at a later date
    atlas=atlasORroiDir;
elseif strcmp(atlasORroiDir(end-5:end),'.nii.gz')
    atlas=niftiRead(atlasORroiDir);
else
    %do nothing
end

outputROIS=[];

%iterate across input requests
for iROIsOut=1:length(multipleRequestedAmalgums)
    %roiDirCase, shouldn't be char at this point if it isn't a folder
    if ischar(atlasORroiDir)
        if isfolder(atlasORroiDir)
            currentRoisRequested=splitlines(compose(strrep(multipleRequestedAmalgums{iROIsOut},' ','\n')));
            %drop the empty cells
            currentRoisRequested=currentRoisRequested(~cellfun(@isempty,(currentRoisRequested)));
            for iROIrequest=1:length(currentRoisRequested)
                currentRequest=currentRoisRequested{iROIrequest};
                %avoiding indexing the end of the string, will just use
                %contains
                if contains(currentRequest,'nii.gz')
                    currentRequestFname=fullfile(atlasORroiDir,currentRequest);
                else %in the case that nii.gz isn't on the end of the name
                    currentRequestFname=fullfile(atlasORroiDir,strcat(currentRequest,'.nii.gz'));
                end
                %now try and load it... maybe
                if isfile (currentRequestFname)
                    currentROI=niftiRead(currentRequestFname);
                else
                    error('file %s does not exist',currentRequestFname)
                end
                
                %now merge them if the input specified to do so
                %kind of elegant, actually
                if iROIrequest==1
                    outputROIS{iROIsOut}=currentROI;
                else
                    outputROIS{iROIsOut}=bsc_mergeROIs(outputROIS{iROIsOut},currentROI);
                end
                outputROIS{iROIsOut}.data=smooth3(outputROIS{iROIsOut}.data,'box',smoothKernel);
                
            end
            %is char but isn't a folder?
        end
        %if the atlas exists and has been loaded
    elseif exist('atlas','var')
        %extract current label requests
        if or(ischar(multipleRequestedAmalgums{iROIsOut}),isstr(multipleRequestedAmalgums{iROIsOut}))
            currentAtlasLabelsRequested=str2num(multipleRequestedAmalgums{iROIsOut});
        elseif isnumeric(multipleRequestedAmalgums{iROIsOut})%if its a number
            currentAtlasLabelsRequested=multipleRequestedAmalgums{iROIsOut};
        else
            multipleRequestedAmalgums{iROIsOut}
            error('Input not recognized')
        end
        % because of how bsc_roiFromAtlasNums handles inputs, we dont need
        % to split lines or do anythign complex, just passing in an integer
        % vector should be fine
        currentROI =dtiRoiFromNiftiObjectSmoothWrapper(atlas,currentAtlasLabelsRequested, smoothKernel);
        
        outputROIS{iROIsOut}=currentROI;
        
    else %if you don't recognize the input
        fprintf('\n %s',atlasORroiDir)
        error('Input not recognized')
    end
end

end