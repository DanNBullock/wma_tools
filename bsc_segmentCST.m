function [classificationOut] =bsc_segmentCST(wbfg,atlas,varargin)
% This function automatedly segments the cortico spinal tract
% from a given whole brain fiber group using the subject's 2009 DK
% freesurfer parcellation. 

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory
% -varargin: priors from previous steps

% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University

%% parameter note & initialization

%create left/right lables
sideLabel={'left','right'};

%[categoryPrior] =bsc_streamlineCategoryPriors_v4(wbfg, fsDir,2)
[inflatedAtlas] =bsc_inflateLabels(atlas,2);

categoryPrior=varargin{1};
%categoryPrior=categoryPrior{1};

%[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

%atlas=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');
cbROINums=[7 8 46 47];
lentiLut=[12 13 51 52];
palLut=[13;52];
wmLut=[2;41];
thalLut=[10;49];

%iterates through left and right sides
for leftright= [1,2]
    %sidenum is basically a way of switching  between the left and right
    %hemispheres of the brain in accordance with freesurfer's ROI
    %numbering scheme. left = 1, right = 2
    sidenum=10000+leftright*1000;
    
     thalAntPlane=bsc_planeFromROI_v2(thalLut(leftright), 'anterior',atlas);
    
    CBRoi=bsc_roiFromAtlasNums(inflatedAtlas,cbROINums ,1);

    MotorROI=bsc_roiFromAtlasNums(inflatedAtlas,[168 128 146 129 170 103 ]+sidenum ,1);
    SpineROI=bsc_roiFromAtlasNums(inflatedAtlas,16 ,1);
    if leftright==1
        otherWM=bsc_roiFromAtlasNums(atlas,41 ,1);
    else
        otherWM=bsc_roiFromAtlasNums(atlas,2 ,1);
    end
    
    [cst, cstBool]=wma_SegmentFascicleFromConnectome(wbfg, [{MotorROI} {SpineROI} {otherWM} {CBRoi} {thalAntPlane}], {'endpoints','and', 'not','not','not'}, 'dud');
    
    nonSubCortNames=cellfun(@isempty,strfind(categoryPrior.names,'subcortical'));
    nonCebNames=cellfun(@isempty,strfind(categoryPrior.names,'cerebellum'));
    nonSubCortBool=ismember(categoryPrior.index,find(nonSubCortNames));
    nonCebBoolBool=ismember(categoryPrior.index,find(nonCebNames));
    

    
    classificationOut=bsc_concatClassificationCriteria(classificationOut,strcat(sideLabel{leftright},'CST'),nonSubCortBool,nonCebBoolBool,cstBool);

end
end
    
    