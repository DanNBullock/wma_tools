function [ROInums] =bsc_atlasROINumsFromCoords_v2(atlasNifti,coords,space)
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
% % (C) Daniel Bullock 2017 Bloomington, Indiana
%% begin code

% read in the appropriate aseg niftifile
if isstring(atlasNifti)
    atlasNifti=niftiRead(atlasNifti);
else
    %do nothing
end


% create concatonated vector and library if structure is input
if isstruct(coords)
    chkSize=size(coords{1});
    
    if ~chkSize(1)==3
        coords{1}=rot90(coords{1});
    else
        %do nothing
    end
    beginCat=coords{1};
    chkSize=size(coords{1});
    coordLibrary=ones(1,chkSize(2));
    for iCoordSets=2:length(coords)
        chkSize=size(coords{iCoordSets});
        if ~chkSize(1)==3
            coords{iCoordSets}=rot90(coords{iCoordSets});
        else
        end
        chkSize=size(coords{iCoordSets});
        newCoordLibrary(1:chkSize(2)=iCoordSets;
        coordLibrary=horzcat(coordLibrary,newCoordLibrary);
        beginCat=horzcat(beginCat,coords{iCoordSets});
    end
    inputCoordsStore=coords;
    coords=beginCat;
else
    chkSize=size(coords);
    
    if ~chkSize(1)==3
        coords=rot90(coords1);
    else
        %do nothing
    end
end


switch space
    case 'acpc'
        imgCoords  = floor(mrAnatXformCoords(atlasNifti.qto_ijk, coords))';
        chkSize=size(coords);
        for iCoords=1:chkSize(2)
            ROIIND=atlasNifti.data(imgCoords(1),imgCoords(2),imgCoords(3));
            ROILibrary(iCoords)=ROIIND;
        end
    case 'img'
        imgCoords=floor(coords);
        chkSize=size(coords);
        for iCoords=1:chkSize(2)
            ROIIND=atlasNifti.data(imgCoords(1),imgCoords(2),imgCoords(3));
            ROILibrary(iCoords)=ROIIND;
        end
    case 'mni'
        %bront's ants stuff, then my stuff
        %csv write, transform, csv read
        %actually no, we aren't going to do this here.
end

if exist('coordLibrary','var')
    for iCoordSets=1:length(inputCoordsStore)
        ROInums{iCoordSets}=ROILibrary(iCoordSets==coordLibrary);
    end
else
    ROInums=ROILibrary;
end
    
end
