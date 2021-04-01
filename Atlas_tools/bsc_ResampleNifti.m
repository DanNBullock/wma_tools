function [outT1] = bsc_ResampleNifti(t1File, outMm, bb, xform)
%[t1New] = bsc_mrAnatResampleT1(t1File, outT1FileName, outMm, bb, xform)
%
%  4/26/2020 Note:  This is a parred down version of mrAnatResampleT1
%  https://github.com/vistalab/vistasoft/blob/master/mrAnatomy/VolumeUtilities/mrAnatResampleT1.m
%  That neither loads the input nifti nor save it down (done to expidite
%  processing and to faciliate use as pipeline intermediary).
%
% Resample a t1 to a specified resolution
%
% [t1New, xformNew] = mrAnatResampleT1(t1File, outMm)
%
% Cropped out of mrGrayResampleNiftiClass, as we sometimes want to resample
% a t1 NIFTI without also resampling class files.
%
% Example:
%   nipath  = 't1-1.nii.gz';
%   outMm   = [1 1 1];
%   mrAnatResampleT1(nipath, outpath, outMm);
%
%
% 4/2009: JW
%
%  4/20/26 Dan Bullock

disp('resampling the input nifti...');

% Get the t1
if ischar(t1File)
    t1 = niftiRead(t1File);
else
    t1=t1File;%do nothing
end

% Get the xform from the nifti struct
if notDefined('xform'), xform = t1.qto_xyz; end

if notDefined('bb')
    % Create a bounding box in image space
    bb = [1 1 1; t1.dim(1:3)+1];

    % Convert the bounding box to mm space
    bb = floor(mrAnatXformCoords(xform,bb));
end

% Reslice
[t1New,xformNew] = mrAnatResliceSpm(double(t1.data),inv(xform),bb,outMm,[0 0 0 7 7 7],0);

% Convert to single
t1New = single(t1New);

%solution from vistasoft
outT1=niftiCreate('data',t1New,'qto_xyz',xformNew);

return