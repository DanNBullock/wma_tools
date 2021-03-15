function [tckOut, classification] = wma_TCKandClassFromOldWMC(oldFGclassified)
%DESCRIPTION:
% This function converts between old-style versions of fg_classified to an
% output tck + classification structure
%
% INPUTS:
% -oldFGclassified: an old fg classified structure. A 1 x N structure, for N white matter structures, with
% "names" and "fibers" fields, wherein "names" are just strings and
% "fibers" are is a cell vector containing 3 x node coordinates for
% streamlines.
%
%
% OUTPUTS:
% -oldFGclassified: the combined tck structure for all input fiber groups
%
% -classification:  a standard, simplified wmc structure.
%
%  (C) Daniel Bullock 2021 Bloomington
%% Begin Code

tckOut=dtiFgArrayToFiberGroup(oldFGclassified, 'allWMStrucs');
 
classification.names={tckOut.subgroupNames.subgroupName};
classification.index=tckOut.subgroup;
%
end