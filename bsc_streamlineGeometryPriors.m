function [asymPrior, effPrior] =bsc_streamlineGeometryPriors(wbfg)
%
% [asymPrior effPrior] =bsc_streamlineGeometryPriors(wbfg)
%
%creates a set of clsasification prior structures indicating the asymetry
%and efficiency ratios of streamlines
%
% Inputs:
% -wbfg: a whole brain fiber group structure
%
% Outputs:
%  asymPrior:  standardly constructed classification structure with tract
%  names being bins of asymetry ratios, as binned by the incriment
%  variable.
%
%  effPrior:  standardly constructed classification structure with tract
%  names being bins of efficiency ratios, as binned by the incriment
%  variable.
%
% (C) Daniel Bullock, 2019, Indiana University

%% parameter note & initialization

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

 effPrior=classificationOut;
asymPrior=classificationOut;

incriment=.05;

[costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ]=ConnectomeTestQ_v2(wbfg);

for iRatios =incriment:incriment:1
    asymPrior=bsc_concatClassificationCriteria(asymPrior,strcat(num2str(iRatios-incriment),'_to_',num2str(iRatios),'_asym'),AsymRat>(iRatios-incriment),AsymRat<iRatios);
    effPrior=bsc_concatClassificationCriteria(effPrior,strcat(num2str(iRatios-incriment),'_to_',num2str(iRatios),'_eff'),efficiencyRat>(iRatios-incriment),efficiencyRat<iRatios);
end
end