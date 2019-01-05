function wma_multiSliceAmalgumROIatCoords_BL()
% [ROIs] =wma_multiSliceAmalgumROIatCoords_BL(atlasNifti,coords,space)
%
%  Purpose:  given a series of coordinates, find the atlas ROIs that
%  correspond to them, merge them in to one omnibus ROI, and then
%  sequentially slice it into regions defined by the coordinates.
%  Autodetects the orientation of the corrdinate sequence.
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
%  ROIs: an N +1 long cell structure with ROIs corresponding to
%
% % (C) Daniel Bullock 2018 Bloomington, Indiana
%% Begin code

if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
end


config = loadjson('config.json');
coordSets=splitlines(config.coords);
interpolateFlag=config.interpolateFlag;
%outRoiTotalNumlength=length(coordSets);
for Ilines=1:length(coordSets)
    coordCount(Ilines)=length(strfind(coordSets{Ilines},','))+1;
end
coordsIn=dlmread('sub_wrp_coords.csv',',',1,0);
coordsIn=coordsIn(1:end,1:3);

atlas=fullfile(config.atlas,'parc.nii.gz');

%% gen ROI
%fprintf('Generating sliced ROIs for the following coordinates: \n %s',num2str(coordsIn));

mkdir('rois/')
roiCount=coordCount+1;
%openkeyfile
system('touch keyfile.csv')
fid=fopen('keyfile.csv','wt');


for Ilines=1:length(coordSets)
    if Ilines==1
        coordIndexes=1:coordCount(1);
    else
        coordIndexes=sum(coordCount(1:Ilines-1))+1:sum(coordCount(1:Ilines-1))+coordCount(Ilines);
    end
    
    
    coordsHold=coordsIn(coordIndexes,:)';
    fprintf('\n Generating sliced ROIs for coordinates:')
    coordsHold
    fprintf('corresponding to these coordinates in mni space')
    str2num(strrep(coordSets{Ilines},',',newline))'
    
    [ROIs] =bsc_sliceAmalgumROIatCoords(atlas,coordsHold,'acpc',interpolateFlag)
    for iROIS=1:roiCount(Ilines)
        
        if Ilines==1
            indexLabel=iROIS;
            fprintf('\n Saving %s as %s',ROIs{iROIS}.name,strcat('ROI',num2str(indexLabel)))
            [~,~]=dtiRoiNiftiFromMat(ROIs{iROIS},atlas,strcat('rois/ROI',num2str(indexLabel),'.nii.gz'),1);
        else
            indexLabel=iROIS+sum(roiCount(1:Ilines-1));
            fprintf('\n Saving %s as %s',ROIs{iROIS}.name,strcat('ROI',num2str(indexLabel)))
            [~,~]=dtiRoiNiftiFromMat(ROIs{iROIS},atlas,strcat('rois/ROI',num2str(indexLabel),'.nii.gz'),1);
        end
  fprintf(fid,'%i,%s \n',indexLabel,ROIs{iROIS}.name);
    end
end
fclose(fid);
end