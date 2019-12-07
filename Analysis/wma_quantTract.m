function [tractStat]= wma_quantTract(fg)
%[tractStat]= wma_quantTract(fg)
%
%  This function computes a number of descriptive quantative features fiber tracts
%
%  Inputs
%
%  fg:  A fiber group
%
%  Outputs
%
%  tractStat:  a strcuture with named fields corresponding to the measured
%              values.  Measured values are specific to the input tract
%
%
% (C) Daniel Bullock, 2017, Indiana University
%% Paramater set
%limit on the size of a tract that this code will try
streamVol=100000;

%% Begin Code

%find fg name
tractStat.name=fg.name;

streamsTotal=fg.fibers;

if length(streamsTotal)==0
    fprintf('\n empty tract for %s',fg.name')
    tractStat.avgAsymRat=0;
    tractStat.stDevAsymRat=0;
    tractStat.avgFullDisp=0;
    tractStat.stDevFullDisp=0;
    tractStat.avgefficiencyRat=0;
    tractStat.stDevefficiencyRat=0;
    tractStat.stream_count=0;
    tractStat.length_total=0;
    tractStat.avg_stream_length=0;
    tractStat.stream_length_stdev=0;
    tractStat.volume=0;
    tractStat.volLengthRatio=0;
    tractStat.lengthCounts=0;
    tractStat.endpointVolume1=0;
    tractStat.avgEndpointCoord1=0;
    tractStat.endpointDensity1=0;
    tractStat.avgEndpointDist1=0;
    tractStat.stDevEndpointDist1=0;
    tractStat.endpointVolume2=0;
    tractStat.avgEndpointCoord2=0;
    tractStat.endpointDensity2=0;
    tractStat.avgEndpointDist2=0;
    tractStat.stDevEndpointDist2=0;
    tractStat.midpointVolume=0;
    tractStat.avgMidpointCoord=0;
    tractStat.midpointDensity=0;
    tractStat.avgMidpointDist=0;
    tractStat.stDevMidpointDist=0;
    
else
    
    %compute effeciency measures and place in tractStat structure
    [costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ] = ConnectomeTestQ_v2(fg);
    
    tractStat.avgAsymRat=mean(AsymRat);
    tractStat.stDevAsymRat=std(AsymRat);
    
    tractStat.avgFullDisp=mean(FullDisp);
    tractStat.stDevFullDisp=std(FullDisp);
    
    tractStat.avgefficiencyRat=mean(efficiencyRat);
    tractStat.stDevefficiencyRat=std(efficiencyRat);
    
    
    %compute basic statistics
    tractStat.stream_count=length(streamsTotal);
    tractStat.length_total=sum(streamLengths);
    
    % avgerage and standard deviation for streamline lengths
    tractStat.avg_stream_length=mean(streamLengths);
    tractStat.stream_length_stdev=std(streamLengths);
    
    %compute volume
    volVec=[];
    for istreamlines=1:length(streamsTotal)
        streamLengths(istreamlines)=sum(sqrt(sum(diff(streamsTotal{istreamlines},1,2).^2)));
        volVec=horzcat(volVec,streamsTotal{istreamlines});
        %prevent memory usage from getting too extreme
        if rem(istreamlines,5000)==0
            volVec=   unique(round(volVec'),'rows')';
        end
    end
    
    %finish volume computation
    tractStat.volume=length(unique(round(volVec'),'rows'));
    %kind of like density
    tractStat.volLengthRatio=tractStat.length_total/tractStat.volume;
    
    %get histcount data for plotting lenght histogram
    [counts,edges]=histcounts(streamLengths,'BinWidth',1,'BinLimits',[1,300]);
    
    tractStat.lengthCounts=counts;
    
    
    %% Point CloudStats
    %will only try and do this if the number of streamlines is managable
    if streamVol>tractStat.volume
        
        %reorient fibers so endpoint clouds make sense.
        fg=bsc_reorientFiber(fg);
        
        %cluster the endpoints
        [RASout, LPIout, RASoutEndpoint, LPIoutEndpoint] = endpointClusterProto(fg);
        
        %get the midpoints
        midpoints=[];
        for iFibers=1:length(streamsTotal)
            fiberNodeNum=round(length(streamsTotal{iFibers})/2);
            curStreamline=streamsTotal{iFibers};
            midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
        end
        
        %compute vol stats for RAS endpoint cloud
        tractStat.endpointVolume1=length(unique(round(RASout),'rows'));
        tractStat.avgEndpointCoord1=mean(RASout,1);
        tractStat.endpointDensity1=tractStat.stream_count/tractStat.endpointVolume1;
        
        %compute endpoint distances from RAS cloud centroid
        for istreamlines=1:length(streamsTotal)
            endpointDists1(istreamlines)=sum(sqrt(sum(diff([tractStat.avgEndpointCoord1',RASout(istreamlines,:)'],1,2).^2)));
        end
        
        %compute average RAS endpoint distance from centroid
        tractStat.avgEndpointDist1=mean(endpointDists1);
        tractStat.stDevEndpointDist1=std(endpointDists1);
        
        
        
        %compute vol stats for LPI endpoint cloud
        tractStat.endpointVolume2=length(unique(round(LPIout),'rows'));
        tractStat.avgEndpointCoord2=mean(LPIout,1);
        tractStat.endpointDensity2=tractStat.stream_count/tractStat.endpointVolume2;
        
        %compute endpoint distances from RAS cloud centroid
        for istreamlines=1:length(streamsTotal)
            endpointDists2(istreamlines)=sum(sqrt(sum(diff([tractStat.avgEndpointCoord2',LPIout(istreamlines,:)'],1,2).^2)));
        end
        
        %compute average RAS endpoint distance from centroid
        tractStat.avgEndpointDist2=mean(endpointDists2);
        tractStat.stDevEndpointDist2=std(endpointDists2);
        
        
        
        %compute vol stats for midpoint cloud
        tractStat.midpointVolume=length(unique(round(midpoints),'rows'));
        tractStat.avgMidpointCoord=mean(midpoints,1);
        tractStat.midpointDensity=tractStat.stream_count/tractStat.midpointVolume;
        
        %compute endpoint distances from RAS cloud centroid
        for istreamlines=1:length(streamsTotal)
            midpointDists(istreamlines)=sum(sqrt(sum(diff([tractStat.avgMidpointCoord',midpoints(istreamlines,:)'],1,2).^2)));
        end
        
        %compute average RAS endpoint distance from centroid
        tractStat.avgMidpointDist=mean(midpointDists);
        tractStat.stDevMidpointDist=std(midpointDists);
        
    else
        warning('\n Endpoint and midpoint stats not run due to overly large tract')
    end
end

end