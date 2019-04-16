function  [LPINii, RASNii] = wma_endpointMapsDecay_v5(ind, classification, fe, wbFG,fsDir, kernelSize, filter)
% function  wma_endpointMapsDecay(fg, fsDir, thresholdinmm, decayFunc)
%
% OVERVIEW: generate .nifti files for (both) of the tract endpoint density
% mappings for a given fibergroup and save them down in the designated
% output directory.
%
% INPUTS:
% -fg: fiber group structure, IN AN IMAGE SPACE COORDINATE SCHEME, as a transform
% to an img space coordinate scheme is later performed
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -saveDir: The designared output directory
% -kernelSize: the smoothKernel to be used in smooth3
% -filter: smooth3 kernel to be used
%          ->box
%          ->gaussian
%
%  Maybe impliment more options via :
%  https://www.mathworks.com/help/stats/kernel-distribution.html
%
% REQUIREMENT: Parallel computing toolbox and MATLAB newer than 2013b.
%
% OUTPUTS:
% -LPINii, RASNii: Endpoint density file for two distinct endpoint clouds of the tract.
% In this file, we simply sum the values across streamline endpoints. LPI =
% left posterior inferior endpoint group, RAS = Right anterior superior
% endpoint bgroup
%
% (C) Daniel Bullock and Hiromasa Takemura, 2016, CiNet HHS
%  Overhaul by DNB 10/2017

%% Parameter settings & preliminaries
tic
% Define a distance threshold. Default is 3 mm.
if notDefined('kernelSize'), kernelSize=3;end

% Default is to compute sum of all endpoints within threshold distance
if notDefined('filter'), filter='box';end

%generate the gm mask if it doesn't exist
if ~exist(strcat(fsDir,'/gm_mask.nii.gz'))
    %apaprently necessary for matlab?
    spaceChar={' '};
    %why are these structures?
    string1=strcat('mri_binarize --i',spaceChar,fsDir,'/mri/aparc+aseg.mgz --gm --o',spaceChar, fsDir, '/gm_mask.mgz');
    string2=strcat('mri_convert',spaceChar,fsDir,'/gm_mask.mgz',spaceChar, fsDir, '/gm_mask.nii.gz');
    [status1 result] = system(string1{1});
    [status2 result] = system(string2{1});
    if or(status1,status2)
        fprintf('/n Error generating gm nifti file.  There may be a problem finding the aaparc+aseg.mgz or the gm_mask.mgz file.')
        keyboard
    end
end

%% Load files

% Load Gray matter mask; previously this was either a left or right mask,
% but now we make a whole brain grey matter mask.  As such we can probably
% handle cross-hemispheric fibers.

%% Estimate the orientation of the tract and get endpoint clouds
%REPLACE WITH THIS FUNCTION
if notDefined('wbFG')
[wbFG, fe] = bsc_LoadAndParseFiberStructure(fe);
end


fg=wbFG;
fg.fibers=wbFG.fibers(classification.index==ind);

% fgIMG=fg;
% fgIMG.fibers=fe.fg.fibers(ind);


[RASout, LPIout, RASoutEndpoint, LPIoutEndpoint] = endpointClusterProto(fg);
% [RASoutIMG, LPIoutIMG, RASoutEndpoint, LPIoutEndpoint] = endpointClusterProto(fgIMG);



graynii = niftiRead(strcat(fsDir,'/gm_mask.nii.gz'));
% 
% grayROI = wma_dtiImportRoiFromNifti(graynii);
% 
% figure
% bsc_plotROIEndpointsOnFG(fg,grayROI)
% figure
% bsc_plotROIEndpointsOnFG(fgIMG,grayROI)


[RASout, ~, ~] = mrAnatXformCoords(graynii.qto_ijk, RASout);
[LPIout, ~, ~] = mrAnatXformCoords(graynii.qto_ijk, LPIout);



%% Set output file names
fname{1} = strcat(strrep(classification.names{ind},' ','_'),'_',filter,'_',num2str(kernelSize),'mm_RAS_FiberEndpoint.nii.gz'); % Tract endpoint density without normalization
fname{2} = strcat(strrep(classification.names{ind},' ','_'),'_',filter,'_',num2str(kernelSize),'mm_LPI_FiberEndpoint.nii.gz');

%% Extract fiber endpoint coordinates



RASNiiData = zeros(graynii.dim(1:3));
LPINiiData = zeros(graynii.dim(1:3));


RASindex=round(RASout)+1;
LPIindex=round(LPIout)+1;


for istreams=1:length(fg.fibers)
    RASNiiData(RASindex(istreams,1),RASindex(istreams,2),RASindex(istreams,3))=RASNiiData(RASindex(istreams,1),RASindex(istreams,2),RASindex(istreams,3))+1;
    LPINiiData(LPIindex(istreams,1),LPIindex(istreams,2),LPIindex(istreams,3))=LPINiiData(LPIindex(istreams,1),LPIindex(istreams,2),LPIindex(istreams,3))+1;
end



%% Apply Smooth

RASNiiData=smooth3(RASNiiData,filter,kernelSize);
LPINiiData=smooth3(LPINiiData,filter,kernelSize);

%% Apply Mask
%not going to do it for now and see what happens
% RASNiiData(~graynii.data)=0;
% LPINiiData(~graynii.data)=0;


%% generate blank nii structures for the endpoint mappings
%naming convention reversed, but it has to be.
% 
% RASNii = niftiCreate('data', RASNiiData, 'fname', fname{1}, ...
%                    'qto_xyz', fe.life.xform.img2acpc, ...
%                    'qto_ijk', fe.life.xform.acpc2img);
% LPINii = niftiCreate('data', LPINiiData, 'fname', fname{2}, ...
%                    'qto_xyz', fe.life.xform.img2acpc, ...
%                    'qto_ijk', fe.life.xform.acpc2img);               

RASNii = graynii;
RASNii.fname = fname{1};
RASNii.data = RASNiiData;

% RASROI = wma_dtiImportRoiFromNifti(RASNii);
% figure
% bsc_plotROIEndpointsOnFG(fg,RASROI)
% figure
% bsc_plotROIEndpointsOnFG(fgIMG,RASROI)




LPINii = graynii;
LPINii.fname = fname{2};
LPINii.data = LPINiiData;

computeTime=toc;
fprintf('\n Density mesh generation for %s with %i streamlines complete in %4.2f seconds',classification.names{ind}, length(fg.fibers),computeTime)

%% set names appropriately

end

