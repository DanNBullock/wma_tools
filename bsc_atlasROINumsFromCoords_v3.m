function [ROInums] =bsc_atlasROINumsFromCoords_v3(atlasNifti,coords,space,interpolateBool)
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

ROILibrary=zeros(length(coords),1);


% eventually use imgspace coords to get appropriate label
switch space
    case 'acpc'

        
        imgCoords  = floor(mrAnatXformCoords(atlasNifti.qto_ijk, coords))';
        chkSize=size(coords);
        for iCoords=1:chkSize(2)
            candidateLabels=unique(atlasNifti.data(imgCoords(1,iCoords)-1:imgCoords(1,iCoords)+1,imgCoords(2,iCoords)-1:imgCoords(2,iCoords)+1,imgCoords(3,iCoords)-1:imgCoords(3,iCoords)+1));
            
            for iCandidates=1:length(candidateLabels)
                labelCount(iCandidates)=sum(sum(sum(atlasNifti.data(imgCoords(1,iCoords)-1:imgCoords(1,iCoords)+1,imgCoords(2,iCoords)-1:imgCoords(2,iCoords)+1,imgCoords(3,iCoords)-1:imgCoords(3,iCoords)+1)==candidateLabels(iCandidates))));
            end
            
            tiebreakCount=1;
            if length(candidateLabels(max(labelCount)==labelCount))>1
                
                %tiebreaker
                
                while length(candidateLabels(max(labelCount)==labelCount))>1
                    %will crap out of this happens on an edge tiebreakCount
                    candidateLabels=unique(atlasNifti.data(imgCoords(1,iCoords)-tiebreakCount:imgCoords(1,iCoords)+tiebreakCount,imgCoords(2,iCoords)-tiebreakCount:imgCoords(2,iCoords)+tiebreakCount,imgCoords(3,iCoords)-tiebreakCount:imgCoords(3,iCoords)+tiebreakCount));

                    clear labelCount
                    for iCandidates=1:length(candidateLabels)
                         labelCount(iCandidates)=sum(sum(sum(atlasNifti.data(imgCoords(1,iCoords)-tiebreakCount:imgCoords(1,iCoords)+tiebreakCount,imgCoords(2,iCoords)-1:imgCoords(2,iCoords)+tiebreakCount,imgCoords(3,iCoords)-tiebreakCount:imgCoords(3,iCoords)+tiebreakCount)==candidateLabels(iCandidates))));
             end
                    tiebreakCount=tiebreakCount+1;
                end
                
                %fprintf('\n tiebreaker resolved after %i iterations for %i voxels',tiebreakCount-1,length(find(atlasIterCur.data(voxelMask))))
                thisLabel=candidateLabels(max(labelCount)==labelCount);
                % thisLabel=0;
                ROILibrary(iCoords)=thisLabel;
  
                clear labelCount
            else
                thisLabel=candidateLabels(max(labelCount)==labelCount);
                ROILibrary(iCoords)=thisLabel;
                clear labelCount
                
            end
            labelCount=[];
            
        end
    case 'img'
        imgCoords=floor(coords);
        chkSize=size(coords);
        for iCoords=1:chkSize(2)
            ROIIND=atlasNifti.data(imgCoords(1,iCoords),imgCoords(2,iCoords),imgCoords(3,iCoords));
            ROILibrary(iCoords)=ROIIND;
            imgSpaceCoord=num2str([imgCoords(1,iCoords),imgCoords(2,iCoords),imgCoords(3,iCoords)]);
            % fprintf('\n %s corresponds to number %i in %s',imgSpaceCoord,ROIIND,path2);
        end
end

%set variable
ROInums=ROILibrary;
end
