function [pool, poolDir]=startUniqueParpool(nclust)
% function pool=startUniqueParpool(nclust)
%% start parallel pool

poolobj = gcp('nocreate');
if sum(size(poolobj))==0
    tmp_folder = '~/tmp/';
    
    
    
    display(['Opening parallel pool with ', num2str(nclust), ' cores...']);
    
    % create parallel cluster object
    
    c = parcluster;
    
    % set number of cores from arguments
    
    c.NumWorkers = nclust;
    
    % set temporary cache directory
    
    t = tempname(tmp_folder);
    
    % make cache dir
    
    OK = mkdir(t);
    
    % check and set cachedir location
    
    if OK
        
        % set local storage for parpool
        
        c.JobStorageLocation = t;
        poolDir=t ;
    end
    
    % catch cache dir path for output
    
    %files.cache = t;
    
    % start parpool - close parpool at end of fxn
    
    pool = parpool(c, nclust, 'IdleTimeout', 120);
else
    fprintf('\n Parpool already exists')
    pool=gcp('nocreate');
    poolDir=[];
end
end