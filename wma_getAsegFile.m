function asegFile = wma_getAsegFile(fsDir , asegOption)
%
% asegFile = wma_getAsegFile(fsDir , asegOption)
%
% Description:  This funcion loads the aseg file for the given path and
% option.  If the aseg file doesn't exist it creates the .mat variant.
%
% Inputs: 
%   -fsDir: path to a particular subject's freesurfer directory
%
%   -asegOption: either 'orig' or '2009'.  
%
% Outputs:
% -asegFile: the aseg .mat structure
%
%  (C) Daniel Bullock 2017 Bloomington

switch asegOption
    case '2009'
        asegName='aparc.a2009s+aseg';
    case 'orig'
        asegName='aparc+aseg';
end  

if ~exist(strcat(fsDir,'/mri/',asegName,'.nii.gz'),'file')
    warning('%s not found, attempting to create file',strcat(fsDir,'/mri/',asegName,'.nii.gz')) 
    %apaprently necessary for matlab?
    spaceChar={' '};
    cmndString=strcat('mri_convert',spaceChar,fsDir,'/mri/',asegName,'.mgz',spaceChar, fsDir, '/mri/',asegName,'.nii.gz');

    [status result] = system(cmndString{1},'-echo');
    if status~=0
        error('/n Error generating aseg nifti file.  There may be a problem finding the .mgz file.  Ensure mri_convert is loaded. Output: %s',result)
    end
end

%reads in label data
asegFile=niftiRead(strcat(fsDir,'/mri/',asegName,'.nii.gz'));

end