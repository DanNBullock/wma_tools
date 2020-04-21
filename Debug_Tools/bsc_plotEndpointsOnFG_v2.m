function bsc_plotEndpointsOnFG_v2(fg)
% bsc_plotEndpointsOnFG_v2(fg)
%
% Purpose:  this function is designed to aid segmentation efforts. it plots
% a collection of streamlines (entered as an FG) along with their endpoints
% and midpoints.  Using this information a user can hone their segmenation.
%
%  V2 adds 50 unit long colored lines for axis reference, colors consistent
%  with standard DWI color scheme
%
% Inputs:
%
% -fg: an fg structure corresponding to a tract as from vistasoft
%
% Outputs:
%
%  [none] - function is made for interactive use and creates a plot, not
%  for general plotting / figure generation use
%
% %  (C) Daniel Bullock 2018 Bloomington
%
%% Begin code

if ~isempty( fg.fibers)
    %We'll be adding color coordinated markers to end, so we need to group
    %the streamlines appropriately
    [cluster1, cluster2, ~, ~]=endpointClusterProto(fg);
    
    %make the reference cross
    figure
    hold on
    %plot x line
    plot3([50 -50],[0 0],[0 0],'r','LineWidth',1)
    %plot y line
    plot3([0 0],[50 -50],[0 0],'g','LineWidth',1)
    %plot z line
    plot3([0 0],[0 0],[50 -50],'b','LineWidth',1)
    
    bsc_quickPlot(fg,'m')
    hold on
    
    %get the midpoints, for labeling
    for iFibers=1:length(fg.fibers)
        fiberNodeNum=round(length(fg.fibers{iFibers})/2);
        curStreamline=fg.fibers{iFibers};
        midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    end
    
    %plot midpoints
    midpoints1=scatter3(midpoints(:,1),midpoints(:,2),midpoints(:,3),3,'k','*');
    
    %plot group 1 endpoints
    endpoints1=scatter3(cluster1(:,1),cluster1(:,2),cluster1(:,3),3,'r');
    
    %plot group 2 endpoints
    endpoints2=scatter3(cluster2(:,1),cluster2(:,2),cluster2(:,3),3,'b');
    
    legend([midpoints1 endpoints1 endpoints2],'midpoints','RAS endpoints','LPI endpoints')
else
    %not much else you can do
    fprintf('\n no fibers to plot for %s', fg.name)
end

end