function [olab] = fnDeislandLabels_v2(labels,outfile,maxisleSize,replaceVal)

% voxles from outside of the largest continuous label (bad freesurfer labels)
%
%   Detailed explanation goes here
%
% Brent McPherson, (c) 2019, Indiana University
% Dan Bullock, slight modifications.
%


% read in aligned aparc+aseg .nii (read labels)
if or(ischar(labels),isstr(labels))
    labs = niftiRead(labels);
else
    labs=labels;
end

% create a copy in output
olab = labs;
if ~notDefined('outfile')
    olab.fname = outfile;
end

fprintf('fnDeislandLabels_v2 max:%i relace:%i\n', maxisleSize, replaceVal)

% get list of unique labels (0 is background)
ulab = unique(labs.data(:));
ulab = ulab(ulab > 0);

% extract the data
data = labs.data;
imgdim = size(labs.data);

% preallocate output data labels
odat = zeros(imgdim);

nlost=0;
for ii = 1:size(ulab, 1)
    
    % create 3d image of the label
    tid = int32(data == ulab(ii));
    
    % find islands of label
    cc = bwareaopen(tid,maxisleSize);
   
    % write the output labels from the largest component
    odat(cc)=ulab(ii);
    
    odat(tid&~cc)=replaceVal;
    nlost=nlost+length(find(tid&~cc));
end
fprintf('%i total voxels lost',nlost)
% assign output data
olab.data = odat;

disp('Saving fixed label nifti...');

% write output
if ~notDefined('outfile')
    niftiWrite(olab);
end

end
