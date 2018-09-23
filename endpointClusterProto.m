function [RASout, LPIout, RASoutEndpoint, LPIoutEndpoint] = endpointClusterProto(fg)


for ifibers= 1:length(fg.fibers)
    testStreamline=fg.fibers{ifibers};
    
    endpoint1(ifibers,:)=testStreamline(:,1);
    
    endpoint2(ifibers,:)=testStreamline(:,end);
    
    quarterPoint1Index(ifibers)=round(length(testStreamline)*.25);
    quarterPoint1(ifibers,:)=testStreamline(:,quarterPoint1Index(ifibers));
    
    quarterPoint2Index(ifibers)=round(length(testStreamline)*.75);
    quarterPoint2(ifibers,:)=testStreamline(:,quarterPoint2Index(ifibers));
    
    fourtyIndex(ifibers)=round(length(testStreamline)*.4);
    fourtyPoint(ifibers,:)=testStreamline(:,fourtyIndex(ifibers));
    
    sixtyIndex(ifibers)=round(length(testStreamline)*.6);
    sixtyPoint(ifibers,:)=testStreamline(:,sixtyIndex(ifibers));
    
    midpointIndex(ifibers)=round(length(testStreamline)*.5);
    midPoint(ifibers,:)=testStreamline(:,midpointIndex(ifibers));
    
    if ifibers==1
        cluster1(ifibers,:)=endpoint1(ifibers,:);
        cluster2(ifibers,:)=endpoint2(ifibers,:);
        
        quarterCluster1(ifibers,:)=quarterPoint1(ifibers,:);
        quarterCluster2(ifibers,:)=quarterPoint2(ifibers,:);
        
        sixtyCluster(ifibers,:)=sixtyPoint(ifibers,:);
        fourtyCluster(ifibers,:)=fourtyPoint(ifibers,:);
        
        endpoint1ID(ifibers)=1;
        endpoint2ID(ifibers)=2;
        
        
    elseif ifibers==2
        
        cluster1ToEndpoint1=pdist(vertcat(cluster1(1,:),endpoint1(ifibers,:)));
        cluster2ToEndpoint1=pdist(vertcat(cluster2(1,:),endpoint1(ifibers,:)));
        
        cluster2ToEndpoint2=pdist(vertcat(cluster2(1,:),endpoint2(ifibers,:)));
        cluster1ToEndpoint2=pdist(vertcat(cluster1(1,:),endpoint2(ifibers,:)));
        
        Endpoint1ProxVote=find(min([cluster1ToEndpoint1,cluster2ToEndpoint1])==[cluster1ToEndpoint1,cluster2ToEndpoint1]);
        Endpoint2ProxVote=find(min([cluster1ToEndpoint2,cluster2ToEndpoint2])==[cluster1ToEndpoint2,cluster2ToEndpoint2]);
        
        
        Qcluster1ToEndpoint1=pdist(vertcat(quarterCluster1(1,:),quarterPoint1(ifibers,:)));
        Qcluster2ToEndpoint1=pdist(vertcat(quarterCluster2(1,:),quarterPoint1(ifibers,:)));
        
        Qcluster2ToEndpoint2=pdist(vertcat(quarterCluster2(1,:),quarterPoint2(ifibers,:)));
        Qcluster1ToEndpoint2=pdist(vertcat(quarterCluster1(1,:),quarterPoint2(ifibers,:)));
        
        Qpoint1ProxVote=find(min([Qcluster1ToEndpoint1,Qcluster2ToEndpoint1])==[Qcluster1ToEndpoint1,Qcluster2ToEndpoint1]);
        Qpoint2ProxVote=find(min([Qcluster1ToEndpoint2,Qcluster2ToEndpoint2])==[Qcluster1ToEndpoint2,Qcluster2ToEndpoint2]);
        
        %The naming is off here, but whatever
        sixtyClusterToSixtyPoint=pdist(vertcat(sixtyCluster(1,:),sixtyPoint(ifibers,:)));
        fourtyClusterToSixtyPoint=pdist(vertcat(fourtyCluster(1,:),sixtyPoint(ifibers,:)));
        
        fourtyClusterToFourtyPoint=pdist(vertcat(fourtyCluster(1,:),fourtyPoint(ifibers,:)));
        sixtyClusterToFourtyPoint=pdist(vertcat(sixtyCluster(1,:),fourtyPoint(ifibers,:)));
        
        sixtyProxVote=find(min([sixtyClusterToSixtyPoint,fourtyClusterToSixtyPoint])==[sixtyClusterToSixtyPoint,fourtyClusterToSixtyPoint]);
        fourtyProxVote=find(min([sixtyClusterToFourtyPoint,fourtyClusterToFourtyPoint])==[sixtyClusterToFourtyPoint,fourtyClusterToFourtyPoint]);
        
        
        Agree1Chk=Endpoint1ProxVote==Qpoint1ProxVote==fourtyProxVote;
        Agree2Chk=Endpoint2ProxVote==Qpoint2ProxVote==sixtyProxVote;
        
        ExclusiveCheck=Endpoint1ProxVote~=Endpoint2ProxVote;
        
        if Agree1Chk && Agree2Chk && ExclusiveCheck
            
            
            endpoint1ID(ifibers)=Endpoint1ProxVote;
            endpoint2ID(ifibers)=Endpoint2ProxVote;
        else
            
            
            Option1=cluster1ToEndpoint1+cluster2ToEndpoint2;
            Option2=cluster2ToEndpoint1+cluster1ToEndpoint2;
            
            QOption1=Qcluster1ToEndpoint1+Qcluster2ToEndpoint2;
            QOption2=Qcluster2ToEndpoint1+Qcluster1ToEndpoint2;
            
            tenOption1=sixtyClusterToSixtyPoint+fourtyClusterToFourtyPoint;
            tenOption2=fourtyClusterToSixtyPoint+sixtyClusterToFourtyPoint;
            
            TotalError1=Option1+QOption1+tenOption1;
            TotalError2=Option2+QOption2+tenOption2;
            
            optionFlag=find(min([TotalError1,TotalError2])==[TotalError1,TotalError2]);
            
            if optionFlag==1
                endpoint1ID(ifibers)=1;
                endpoint2ID(ifibers)=2;
            else
                endpoint1ID(ifibers)=2;
                endpoint2ID(ifibers)=1;
                
            end
        end
        
        
        if endpoint1ID(ifibers) == 1
            cluster1(ifibers,:)=endpoint1(ifibers,:);
            cluster2(ifibers,:)=endpoint2(ifibers,:);
            
            quarterCluster1(ifibers,:)=quarterPoint1(ifibers,:);
            quarterCluster2(ifibers,:)=quarterPoint2(ifibers,:);
            
            fourtyCluster(ifibers,:)=fourtyPoint(ifibers,:);
            sixtyCluster(ifibers,:)=sixtyPoint(ifibers,:);
            
        else
            cluster2(ifibers,:)=endpoint1(ifibers,:);
            cluster1(ifibers,:)=endpoint2(ifibers,:);
            
            quarterCluster1(ifibers,:)=quarterPoint2(ifibers,:);
            quarterCluster2(ifibers,:)=quarterPoint1(ifibers,:);
            
            fourtyCluster(ifibers,:)=sixtyPoint(ifibers,:);
            sixtyCluster(ifibers,:)=fourtyPoint(ifibers,:);
        end
        
    else
        
        cluster1ToEndpoint1=pdist(vertcat(mean(cluster1,1),endpoint1(ifibers,:)));
        cluster2ToEndpoint1=pdist(vertcat(mean(cluster2,1),endpoint1(ifibers,:)));
        
        cluster2ToEndpoint2=pdist(vertcat(mean(cluster2,1),endpoint2(ifibers,:)));
        cluster1ToEndpoint2=pdist(vertcat(mean(cluster1,1),endpoint2(ifibers,:)));
        
        Endpoint1ProxVote=find(min([cluster1ToEndpoint1,cluster2ToEndpoint1])==[cluster1ToEndpoint1,cluster2ToEndpoint1]);
        Endpoint2ProxVote=find(min([cluster1ToEndpoint2,cluster2ToEndpoint2])==[cluster1ToEndpoint2,cluster2ToEndpoint2]);
        
        Qcluster1ToEndpoint1=pdist(vertcat(mean(quarterCluster1,1),quarterPoint1(ifibers,:)));
        Qcluster2ToEndpoint1=pdist(vertcat(mean(quarterCluster2,1),quarterPoint1(ifibers,:)));
        
        Qcluster2ToEndpoint2=pdist(vertcat(mean(quarterCluster2,1),quarterPoint2(ifibers,:)));
        Qcluster1ToEndpoint2=pdist(vertcat(mean(quarterCluster1,1),quarterPoint2(ifibers,:)));
        
        Qpoint1ProxVote=find(min([Qcluster1ToEndpoint1,Qcluster2ToEndpoint1])==[Qcluster1ToEndpoint1,Qcluster2ToEndpoint1]);
        Qpoint2ProxVote=find(min([Qcluster1ToEndpoint2,Qcluster2ToEndpoint2])==[Qcluster1ToEndpoint2,Qcluster2ToEndpoint2]);
        
        %The naming is off here, but whatever
        sixtyClusterToSixtyPoint=pdist(vertcat(mean(sixtyCluster,1),sixtyPoint(ifibers,:)));
        fourtyClusterToSixtyPoint=pdist(vertcat(mean(fourtyCluster,1),sixtyPoint(ifibers,:)));
        
        fourtyClusterToFourtyPoint=pdist(vertcat(mean(fourtyCluster,1),fourtyPoint(ifibers,:)));
        sixtyClusterToFourtyPoint=pdist(vertcat(mean(sixtyCluster,1),fourtyPoint(ifibers,:)));
        
        sixtyProxVote=find(min([sixtyClusterToSixtyPoint,fourtyClusterToSixtyPoint])==[sixtyClusterToSixtyPoint,fourtyClusterToSixtyPoint]);
        fourtyProxVote=find(min([sixtyClusterToFourtyPoint,fourtyClusterToFourtyPoint])==[sixtyClusterToFourtyPoint,fourtyClusterToFourtyPoint]);
        
        
        Agree1Chk=Endpoint1ProxVote==Qpoint1ProxVote==fourtyProxVote;
        Agree2Chk=Endpoint2ProxVote==Qpoint2ProxVote==sixtyProxVote;
        
        ExclusiveCheck=Endpoint1ProxVote~=Endpoint2ProxVote;
        
        if Agree1Chk && Agree2Chk && ExclusiveCheck
            
            
            endpoint1ID(ifibers)=Endpoint1ProxVote;
            endpoint2ID(ifibers)=Endpoint2ProxVote;
        else
            
            
            Option1=cluster1ToEndpoint1+cluster2ToEndpoint2;
            Option2=cluster2ToEndpoint1+cluster1ToEndpoint2;
            
            QOption1=Qcluster1ToEndpoint1+Qcluster2ToEndpoint2;
            QOption2=Qcluster2ToEndpoint1+Qcluster1ToEndpoint2;
            
            tenOption1=sixtyClusterToSixtyPoint+fourtyClusterToFourtyPoint;
            tenOption2=fourtyClusterToSixtyPoint+sixtyClusterToFourtyPoint;
            
            TotalError1=Option1+QOption1+tenOption1;
            TotalError2=Option2+QOption2+tenOption2;
            
            optionFlag=find(min([TotalError1,TotalError2])==[TotalError1,TotalError2]);
            
            if optionFlag==1
                endpoint1ID(ifibers)=1;
                endpoint2ID(ifibers)=2;
            else
                endpoint1ID(ifibers)=2;
                endpoint2ID(ifibers)=1;
                
            end
            
        end
            if  endpoint1ID(ifibers) == 1
                cluster1(ifibers,:)=endpoint1(ifibers,:);
                cluster2(ifibers,:)=endpoint2(ifibers,:);
                
                quarterCluster1(ifibers,:)=quarterPoint1(ifibers,:);
                quarterCluster2(ifibers,:)=quarterPoint2(ifibers,:);
                
                fourtyCluster(ifibers,:)=fourtyPoint(ifibers,:);
                sixtyCluster(ifibers,:)=sixtyPoint(ifibers,:);
                
            else
                cluster2(ifibers,:)=endpoint1(ifibers,:);
                cluster1(ifibers,:)=endpoint2(ifibers,:);
                
                quarterCluster1(ifibers,:)=quarterPoint2(ifibers,:);
                quarterCluster2(ifibers,:)=quarterPoint1(ifibers,:);
                
                fourtyCluster(ifibers,:)=sixtyPoint(ifibers,:);
                sixtyCluster(ifibers,:)=fourtyPoint(ifibers,:);
                
            end
        end
        
  
        
        
    end

    


cluster1Centroid=mean(cluster1,1);
cluster2Centroid=mean(cluster2,1);

clusterCentroidDistance=abs(cluster1Centroid-cluster2Centroid);

dimOffset=find(max(clusterCentroidDistance)==clusterCentroidDistance);

RASind=find(max([cluster1Centroid(dimOffset) cluster2Centroid(dimOffset)])==[cluster1Centroid(dimOffset) cluster2Centroid(dimOffset)] );

if RASind==1
    RASout=cluster1;
    RASoutEndpoint=endpoint1ID;
    LPIout=cluster2;
    LPIoutEndpoint=endpoint2ID;
else
        RASout=cluster2;
    RASoutEndpoint=endpoint2ID;
    LPIout=cluster1;
    LPIoutEndpoint=endpoint1ID;


end