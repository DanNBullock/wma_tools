function  bsc_quickPlotClassByName(wbfg,classification,className)
%  bsc_quickPlotClassByName(wbfg,classification,className)
%
%  this function creates a quick and dirty plot of the requested tract.  If
%  no class name is entered, or a match isn't found, a UI prompt will
%  appear. INTENDED FOR INTERACTIVE MATLAB DEBUGGING, NOT FOR COMPILED OR
%  NO INTERFACE USE.
%% begin code

%variable prelims

%if nothing is entered, have a ui pop up.
if isempty(className)
    className = uidropdown(fig,'Items',classification.names,...
        'Value',classification.names{1});
    %now check if it actually matches
elseif isempty(find(strcmp(className,classification.names)))
    className = uidropdown(fig,'Items',classification.names,...
        'Value',classification.names{1});
else %do nothing, it was entered properly, and we can assume its there
end
    
%create a boolean vector just for these streamlines
%no need to reinvent the wheel
%currentStreamsBool=  bsc_extractStreamIndByName(classification,className);
[currentSignleClass]=bsc_extractClassification(classification,className);

%use some old code to make a singleton tract object
tractStruc = bsc_makeFGsFromClassification_v5(currentSignleClass, wbfg);

%go ahead and plot
bsc_plotEndpointsOnFG_v2(tractStruc{1})
end

    
   
