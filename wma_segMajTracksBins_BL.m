function wma_segMajTracksBins_BL()
% [classification] =wma_segMajTracks__BL(wbfg, fsDir)
%
% This function automatedly segments the major human white matter tracts
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation.  Brain-life version

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% classificaton: a standardly organized classification structure
% (C) Daniel Bullock, Indiana University
%% Begin Code
if ~isdeployed
    disp('adding paths');
    addpath(genpath('/N/soft/rhel7/spm/8')) %spm needs to be loaded before vistasoft as vistasoft provides anmean that works
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
    addpath(genpath('/N/soft/rhel7/mrtrix/3.0/mrtrix3/matlab'))
end

config = loadjson('config.json');
wbfg=wma_loadTck(config.track);
%wbfg = dtiImportFibersMrtrix(config.track, .5);

fsDir=strcat(pwd,'/freesurfer');

atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

if ~exist(atlasPath,'file')
    fprintf('\n FILE %i NOT FOUND',atlasPath)
end

classificationHold=[];
classificationHold.names=[];
classificationHold.index=zeros(length(wbfg.fibers),1);
tic

wbfgRemainder=rem(length(wbfg.fibers),500000);
binWalls=0:500000:length(wbfg.fibers);
if ~wbfgRemainder==0
    binWalls=horzcat(binWalls,binWalls(end)+wbfgRemainder);
end

fprintf ('\n %i bins to segment',length(binWalls)-1)
%NEW:  SPLTTING THE FG INTO SUBSECTIONS TO REDUCE MEMORY FOOTPRINT.
for iDivisions=1:length(binWalls)-1
    fprintf('\n segmenting streams %i to %i',binWalls(iDivisions)+1,binWalls(iDivisions+1));
    wbfgSubSection=wbfg;
    wbfgSubSection.fibers=wbfg.fibers(binWalls(iDivisions)+1:binWalls(iDivisions+1));
    
    classificationOut=[];
    classificationOut.names=[];
    classificationOut.index=zeros(length(wbfgSubSection.fibers),1);
    
    fprintf('\n creating priors')
    [categoryPrior] =bsc_streamlineCategoryPriors_v6(wbfgSubSection, fsDir,2);
    [~, effPrior] =bsc_streamlineGeometryPriors(wbfgSubSection);
    
    fprintf('\n prior creation complete')
    
    [AntPostclassificationOut] =bsc_segmentAntPostTracts_v2(wbfgSubSection, fsDir,categoryPrior,effPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,AntPostclassificationOut);
    
    [CCclassificationOut] =bsc_segmentCorpusCallosum_v3(wbfgSubSection, fsDir,0,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,CCclassificationOut);
    
    [SubCclassificationOut] =bsc_segmentSubCortical_v2(wbfgSubSection, fsDir,categoryPrior,effPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,SubCclassificationOut);
    
    [AslantclassificationOut] =bsc_segmentAslant(wbfgSubSection, fsDir,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,AslantclassificationOut);
    
    [MDLFclassificationOut] =bsc_segmentMdLF_ILF_v4(wbfgSubSection, fsDir,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,MDLFclassificationOut);
    
    [pArcTPCclassificationOut] = bsc_segpArcTPC(wbfgSubSection, fsDir);
    classificationOut=bsc_reconcileClassifications(classificationOut,pArcTPCclassificationOut);
    
    [opticclassificationOut] =bsc_opticRadiationSeg_V7(wbfgSubSection, fsDir, 0,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,opticclassificationOut);
    
    [cerebellarclassificationOut] =bsc_segmentCerebellarTracts_v2(wbfgSubSection, fsDir,0,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,cerebellarclassificationOut);
    
    [VOFclassificationOut] =bsc_segmentVOF_v2(wbfgSubSection, fsDir,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,VOFclassificationOut);
    
    [CSTclassificationOut] =bsc_segmentCST(wbfgSubSection, fsDir,categoryPrior);
    classificationOut=bsc_reconcileClassifications(classificationOut,CSTclassificationOut);
    
    classificationHold= bsc_spliceClassifications(classificationHold,classificationOut);
end

classification= wma_resortClassificationStruc(classificationHold);
savepath=strcat(pwd,'classification.mat');
which('classification')
save(savepath,'classification');
which('classification')
toc

wma_formatForBrainLife_v2(classification,wbfg);

end

