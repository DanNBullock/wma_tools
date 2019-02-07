function [results]= wma_quantAllWMNorm(feORwbfg,classification)
%function[figHandle, results]= bsc_feAndSegqualityCheck(fe, classification, saveDir)
%
% This function computes a number of statistics for an input fe structure,
% and, if included, a classification structure.  if you include a 
%
% INPUTS: -feORwbfg:  either a whole brain fiber group / tractogram
% path/object or an fe path / object.
%
% -classification: Either the path to structure or the structure itself.
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
% OUTPUTS:
%
% -results: a structure with multiple fields summarizing various properties
%       of your input FE structure and (if included) Seg segmentation.  See
%       code for specifics of all fields and their interpretation
%
% (C) Daniel Bullock, 2017, Indiana University
%% Begin code

%do whole brain stats
[results]= wma_quantWBFG(feORwbfg);

%extract respective structures
[wbFG, fe] = bsc_LoadAndParseFiberStructure(feORwbfg);

%create fg structures for each individual tract
tractStruc = bsc_makeFGsFromClassification_v4(classification, wbFG);

%get postive fibers
if ~isempty(fe)
    posIndexes=find(fe.life.fit.weights>0);
    classificationPos=wma_clearNonvalidClassifications(classification,fe);
    posTractStruc = bsc_makeFGsFromClassification(classificationPos, wbFG);
end


for itracts=1:length(tractStruc)
    %run stats for indivudal tract
    [tractStatWB]= wma_quantTract(tractStruc{itracts});
    results.WBFG.tractStats{itracts}=tractStatWB;
    
    fprintf('\n %s',tractStruc{itracts}.name)
    
    
    %compute and store normalized values
    results.WBFG.tractStats{itracts}.norms.volumeProp=results.WBFG.tractStats{itracts}.volume/results.WBFG.volume;
    results.WBFG.tractStats{itracts}.norms.countProp=results.WBFG.tractStats{itracts}.stream_count/results.WBFG.stream_count;
    results.WBFG.tractStats{itracts}.norms.wireProp=results.WBFG.tractStats{itracts}.length_total/results.WBFG.length_total;
    if results.WBFG.tractStats{itracts}.volume<100000
    results.WBFG.tractStats{itracts}.norms.endpointVolume1Prop=results.WBFG.tractStats{itracts}.endpointVolume1/results.WBFG.volume;
    results.WBFG.tractStats{itracts}.norms.endpointVolume2Prop=results.WBFG.tractStats{itracts}.endpointVolume2/results.WBFG.volume;
    results.WBFG.tractStats{itracts}.norms.midpointVolumeProp=results.WBFG.tractStats{itracts}.midpointVolume/results.WBFG.volume;
    end
    if ~isempty(fe)
         %run stats for indivudal tract
    [posTractStatWB]= wma_quantTract(posTractStruc{itracts});
    results.LiFEstats.tractStats{itracts}=posTractStatWB;
    
    %compute and store normalized values

    results.LiFEstats.tractStats{itracts}.norms.volumeProp=results.WBFG.tractStats{itracts}.volume/results.LiFEstats.posWBFG.volume;
    results.LiFEstats.tractStats{itracts}.norms.countProp=results.WBFG.tractStats{itracts}.stream_count/results.LiFEstats.posWBFG.stream_count;
    results.LiFEstats.tractStats{itracts}.norms.wireProp=results.WBFG.tractStats{itracts}.length_total/results.LiFEstats.posWBFG.length_total;
     if results.WBFG.tractStats{itracts}.volume<100000
    results.LiFEstats.tractStats{itracts}.norms.endpointVolume1Prop=results.WBFG.tractStats{itracts}.endpointVolume1/results.LiFEstats.posWBFG.volume;
    results.LiFEstats.tractStats{itracts}.norms.endpointVolume2Prop=results.WBFG.tractStats{itracts}.endpointVolume2/results.LiFEstats.posWBFG.volume;
    results.LiFEstats.tractStats{itracts}.norms.midpointVolumeProp=results.WBFG.tractStats{itracts}.midpointVolume/results.LiFEstats.posWBFG.volume;
     end
    end
end
    
end