function [ROInums] =bsc_atlasROINumsFromCoords(atlasNifti,coords,space,interpolateBool)
%[ROInums] =bsc_atlasROINumsFromCoords(atlasNifti,coords,space)
%
%  Purpose:  find the label number associated with a particular coordinate
%  in a given atlas.
%
%  INPUTS:
%  -atlasNifti:  path to an atlas nifti.  An object works too.
%
%  -coords:  a 3 by N (where N is the number of coordinates) vector that
%  you would like to know the atlas ROI numbers for.
%
%  -space:  the space that the coordinates are in, either 'img' or 'acpc'.
%
%  OUTPUTS:
%
%  ROInums: an N long integer vector with the ROI label number for each
%  of the coord input coordinates.
%
% % (C) Daniel Bullock 2018 Bloomington, Indiana
%% begin code

% read in the appropriate aseg niftifile
if or(isstring(atlasNifti),ischar(atlasNifti))
    atlasNifti=niftiRead(atlasNifti);
    [~, path2, ~]=fileparts(atlasNifti.fname);
else
    %do nothing
      [~, path2, ~]=fileparts(atlasNifti.fname);
end

% eventually use imgspace coords to get appropriate label
switch space
    case 'acpc'
        imgCoords  = floor(mrAnatXformCoords(atlasNifti.qto_ijk, coords))';
        chkSize=size(coords);
        for iCoords=1:chkSize(2)
            ROIIND=atlasNifti.data(imgCoords(1,iCoords),imgCoords(2,iCoords),imgCoords(3,iCoords));
            ROILibrary(iCoords)=ROIIND;
            thisCoord=num2str([coords(1,iCoords),coords(2,iCoords),coords(3,iCoords)]);
            imgSpaceCoord=num2str([imgCoords(1,iCoords),imgCoords(2,iCoords),imgCoords(3,iCoords)]);
            fprintf('\n %s (%s) corresponds to number %i in %s',thisCoord,imgSpaceCoord,ROIIND,path2);
        end
    case 'img'
        imgCoords=floor(coords);
        chkSize=size(coords);
        for iCoords=1:chkSize(2)
            ROIIND=atlasNifti.data(imgCoords(1,iCoords),imgCoords(2,iCoords),imgCoords(3,iCoords));
            ROILibrary(iCoords)=ROIIND;
            imgSpaceCoord=num2str([imgCoords(1,iCoords),imgCoords(2,iCoords),imgCoords(3,iCoords)]);
             fprintf('\n %s corresponds to number %i in %s',imgSpaceCoord,ROIIND,path2);
        end
end

%set variable
ROInums=ROILibrary;
end
