function [orientedFiber] =bsc_orientFgRAS_LPI(fg)
% [orientedFiber] =bsc_orientFgRAS_LPI(fg)
%
% This function reorients an input fiber group (presumably reprsenting a
% single tract, this will not do anything useful on a whole brain
% connectome), in such a way that the first node of each streamline is
% associated with the RAS (right, anterior, superior) endpoint cluster and
% the last node is associated with the LPI (left, posterior, inferior)
% endpoint cluster.  The orientation of a fiber group (anterior-posterior,
% left-right, superior-inferior) is determined by the primary dimension of
% traversal, i.e. the dimension in which the tract traverses the most
% space.
%
% Inputs:
% -fg:  the input fiber group.  It is highly recomended that you *do not* put
% a whole brain fiber group in as this input.  The output will not be
% meaningful (there is no "primary dimension of traversal" for a whole
% brain fiber group, and it will take a good deal of time).
%
% Outputs:
%  orientedFiber:  a fg oriented in the manner described in the function
%  description, with first nodes associated with the RAS endpoint group and
%  last nodes associated with the LPI endpoint group.
%
% (C) Daniel Bullock,June 6 2020, Indiana University
%% Begin code

%identify the relevant endpoints and their groups
 [~, ~, RASoutEndpoint, LPIoutEndpoint] = endpointClusterProto(fg);
 
 %extract the sttreams so we can index them
 tractStreams=fg.fibers;
 
 %loop across streamlines
 for iStreams=1:length(tractStreams)
     %if the current streamline oriented such that endpoint 1 is part
     %of the RAS group
      if RASoutEndpoint(iStreams)==1
          %do nothing, because this is what we want
      else %otherwise
          %reorient it such that it is oriented in this fashion
        tractStreams{iStreams}=fliplr(tractStreams{iStreams});
      end
 end
 %now that the streamlines have been reoriented, replace the input object
 %from the fibergoup
 fg.fibers=tractStreams;
 %set the output
 orientedFiber=fg;
 end