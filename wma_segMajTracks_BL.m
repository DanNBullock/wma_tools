function [classificationOut] =wma_segMajTracks_BL()
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
config = loadjson('config.json');

wbfg = dtiImportFibersMrtrix(config.track, .5);

fsDir='freesurfer';

classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);
tic
 [categoryPrior] =bsc_streamlineCategoryPriors_v4(wbfg, fsDir,2);
 [asymPrior, effPrior] =bsc_streamlineGeometryPriors(wbfg);

 [AntPostclassificationOut] =bsc_segmentAntPostTracts(wbfg, fsDir,categoryPrior,effPrior);
 classificationOut=bsc_reconcileClassifications(classificationOut,AntPostclassificationOut);
 
 [CCclassificationOut] =bsc_segmentCorpusCallosum_v3(wbfg, fsDir,0,categoryPrior);
 classificationOut=bsc_reconcileClassifications(classificationOut,CCclassificationOut);
 
 [SubCclassificationOut] =bsc_segmentSubCortical(wbfg, fsDir,categoryPrior,effPrior);
 classificationOut=bsc_reconcileClassifications(classificationOut,SubCclassificationOut);
 
 [AslantclassificationOut] =bsc_segmentAslant(wbfg, fsDir,categoryPrior);
 classificationOut=bsc_reconcileClassifications(classificationOut,AslantclassificationOut);
 
 [MDLFclassificationOut] =bsc_segmentMdLF_ILF_v3(wbfg, fsDir);
 classificationOut=bsc_reconcileClassifications(classificationOut,MDLFclassificationOut);
 
 [pArcTPCclassificationOut] = bsc_segpArcTPC(wbfg, fsDir);
 classificationOut=bsc_reconcileClassifications(classificationOut,pArcTPCclassificationOut);
 
 [opticclassificationOut] =bsc_opticRadiationSeg_V6(wbfg, fsDir, 0);
 classificationOut=bsc_reconcileClassifications(classificationOut,opticclassificationOut);
 
  [cerebellarclassificationOut] =bsc_segmentCerebellarTracts(wbfg, fsDir,0,varargin)
  classificationOut=bsc_reconcileClassifications(classificationOut,cerebellarclassificationOut);
 toc

 wma_formatForBrainLife()


end