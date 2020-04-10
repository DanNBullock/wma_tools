function bsc_plotEndpointsOnFG(fg)
% bsc_plotEndpointsOnFG(fg)
%
% Purpose:  this function is designed to aid segmentation efforts. it plots
% a collection of streamlines (entered as an FG) along with their endpoints
% and midpoints.  Using this information a user can hone their segmenation
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
    [cluster1, cluster2, endpoint1ID, endpoint2ID]=endpointClusterProto(fg);
    
    bsc_quickPlot(fg,'g')
    hold on
    
    for iFibers=1:length(fg.fibers)
        fiberNodeNum=round(length(fg.fibers{iFibers})/2);
        curStreamline=fg.fibers{iFibers};
        midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    end
    
    scatter3(midpoints(:,1),midpoints(:,2),midpoints(:,3),3,'m','*')
    
    scatter3(cluster1(:,1),cluster1(:,2),cluster1(:,3),3,'r')
    
    scatter3(cluster2(:,1),cluster2(:,2),cluster2(:,3),3,'b')
else
    fprintf('\n no fibers to plot for %s', fg.name)
    
end