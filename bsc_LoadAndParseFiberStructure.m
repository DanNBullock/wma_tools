function [wbFG, fe] = bsc_LoadAndParseFiberStructure(feORwbfg)
% 
% [wbFG, fe] = bsc_LoadAndParseFiberStructure(feORwbfg)
%
% Purpose:  Takes in an input, be it a string or an object, loads it if
% needed and then extracts the wbfg in acpc space (if an fe is passed).  An
% attempt at making functions more robust.  Unknown how to handle cases in
% which an img coordspace wbFG is passed.
%
% INPUTS
%
% -feORwbfg: either a string or an object, to either a wbFG or an FE
%  structure
%
% OUTPUTS
%
% -wbFG:  a wbFG structure
%
% -fe: an fe Structure
%
% (C) Daniel Bullock, 2017, Indiana University

% checks if passed variable is a path or an object
if ischar(feORwbfg)
    % loads if a path
    if strcmp(feORwbfg(end-3:end),'.tck')
        wbFG =  dtiImportFibersMrtrix(feORwbfg, .5);
        fe=[];
    else
        load(feORwbfg)
        % function is agnostic about naming conventions, but utilizes
        % regularity in structure fields in order to identify variables of
        % interest.
        varList=whos;
        % for each variable in the workspace, check and see if it is either a
        % life structure (if 'life' field is present) or an fg (if 'fibers'
        % field is present.
        for iVars=1:length(varList)
            curVar=eval(varList(iVars).name);
            if isfield(curVar, 'life')
                fprintf('\n fe structure loaded \n')
                % obtains fibers in acpc space
                wbFG = feGet(curVar, 'fibers acpc');
                fprintf('\n streamlines extracted to wbfg in acpc coordspace \n')
                fe=curVar;
            elseif isfield(curVar,'fibers')
                fprintf('\n fg structure loaded.  \n fg structure is in %s coordpace \n', curvar.coordspace)
                wbFG = curVar;
                fe=[];
            end
        end
    end
    %if an object is passed, proceeds accordingly
elseif isstruct(feORwbfg)
    if isfield(feORwbfg, 'life')
        fprintf('\n fe structure loaded')
        % obtains fibers in acpc space
        wbFG = feGet(feORwbfg, 'fibers acpc');
        fe=feORwbfg;
    elseif isfield(feORwbfg,'fibers')
        if isfield(feORwbfg,'coordspace')
        fprintf('\n fg structure loaded.  \n fg structure is in %s coordpace \n', feORwbfg.coordspace)
        end
        wbFG = feORwbfg;
        fe=[];
    end
else
    % errors if neither a structure nor a string is passed.  Dont know how
    % this would happen in typical use.
    error('\n input is neither a structure nor a string.\n')
end

end
