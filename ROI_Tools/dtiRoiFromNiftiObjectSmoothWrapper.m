function roiOut= dtiRoiFromNiftiObjectSmoothWrapper(nifti,maskValue,smoothKernel)
%roiOut= dtiRoiFromNiftiObjectSmoothWrapper(nifti,maskValue,smoothKernel)
%
% A more reliable replacement for bsc_roiFromAtlasNums
%
% Inputs
%
% nifti: a nifti object featuring integer values in the .data field
%
% maskValue:  those labels that you wish to extract as a binarized roi
%
% smooth kernel:  the smooth kernel you'd like to apply to the output mask.
% This can be used to inflate the anatomy represented by the output
% roi/mask.
%
% Outputs
% roiOut: an roi mask coresponding to the requested labels from input
% variable maskValue
% 
% 4/22/20
% sheepishly adapted from vistasoft by Dan Bullock
%
%% begin code

%to make names readable
maskName=strrep(num2str(maskValue),'  ','_');
%use it once to make an uninflated nifti
unsmoothedNifti = dtiRoiFromNiftiObject(nifti,maskValue,maskName,'nifti',true,false);

%establish object
smoothedNifti=unsmoothedNifti;

%smooth data, and accept all non zero entries
smoothedNifti.data=~(smooth3(smoothedNifti.data,'box',smoothKernel))==0;

%compute percentage increase
increasePct=((length(find(smoothedNifti.data))/length(find(unsmoothedNifti.data)))-1)*100;
fprintf('ROI size increased by %4.2f percent\n',increasePct)

%now use the same function to make it into a .mat roi
roiOut = dtiRoiFromNiftiObject(smoothedNifti,1,maskName,'mat',true,false);
 
end