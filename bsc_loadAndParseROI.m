function [ROI] = bsc_loadAndParseROI(ROIorNifti)
% [ROI] = bsc_loadAndParseROI(ROIorNifti)
%
% Purpose:  Takes in an input, be it a string or an object, loads it if
% needed and then perform appropriate operations to obtain a vistasoft-relevant roi. 
%
% INPUTS
%
% -ROIorNifti: either a string or an object, to either a wbFG or an FE
%  structure
%
% OUTPUTS
%
% -ROI:  an roi (as if from dtiNewRoi)
%
%
% (C) Daniel Bullock, 2018, Indiana University
%
%% Begin code

if  isstring (ROIorNifti)
    %decomposes file path
    [fpath,fname,fext]=fileparts(ROIorNifti);
    if or(strcmp(fext,'.gz'),strcmp(fext,'.nii'))
        niftiIn=niftiRead(ROIorNifti);
        ROI=dtiRoiFromNifti(niftiIn,0,'ROI','mat',true,false);
    else
        %could cause a problem due to how old vistasoft tends to throw up a
        %UI prompt in lack of input cases.
        ROI=dtiReadRoi(ROIorNifti);
    end
elseif isstruct(ROIorNifti)
    if ~isfield(ROIorNifti,'coords')
        ROI=dtiRoiFromNifti(ROIorNifti,0,'ROI','mat',true,false);
    elseif isfield(ROIorNifti,'coords')
        ROI=ROIorNifti;
    end
end