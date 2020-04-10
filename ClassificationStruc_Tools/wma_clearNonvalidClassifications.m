function classification=wma_clearNonvalidClassifications(classification,fe)
%
% FiberIndexes= wma_clearNonvalidClassifications(FiberIndexes,wbFG)
%
% Performs outlier removal given input parameters. Defauts to 4 and 4 if
% nothing put in.
%
% INPUTS: 
%
% -classification: Either the path to structure or the structure itself.
%  The strucure has a field "names" with (N) names of the tracts classified
%  while the field "indexes" has a j long vector (where  j = the nubmer of
%  streamlines in wbFG (i.e. length(wbFG.fibers)).  This j long vector has
%  a 0 for to indicate a streamline has gone unclassified, or a number 1:N
%  indicatate that the streamline has been classified as a member of tract
%  (N).
%
% -wbFG: Either the path to the fe, or the structure itself
%
% considering how short this is, it probaby isn't necessary, but whatever.
% (C) Daniel Bullock 2017 Bloomington

%% preliminaries
if ischar(classification)
    load(classification);
end

if ischar(fe)
    load(fe);
end

% throw error early if improper classification indexing is found.  If the
% number of streamlines classified in classification.index isn't equal to
% the number of streamlines in the fe.fg.fibers structure then we can't be
% sure that the classificatinon for streamline X (of classification.index)
% applies streamline X (of fe.fg.fibers)
if length(classification.index)~=length(fe.fg.fibers)
    warning('\N classification - fe structure mismatch.  Check streamline totals')
end

%% set nonvalidated classifications to 0
validIndices=find(fe.life.fit.weights>0);
fprintf('\n %i positively weighted fibers found in fe structure', length(validIndices));

blankIndex=1:length(fe.life.fit.weights);
invalidIndices=setdiff(blankIndex,validIndices);
fprintf ('\n %i pre-life streamlines classified \n', sum (~classification.index==0))
classification.index(invalidIndices)=0;
fprintf ('\n %i post-life streamlines classified \n', sum (~classification.index==0))
end
