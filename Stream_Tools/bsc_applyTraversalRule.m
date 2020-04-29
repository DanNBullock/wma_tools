function booleanOut=bsc_applyTraversalRule(wbfg,criteriaSequence,relativePositioning)




wbfg=bsc_reorientFiber(wbfg);
allStreams=wbfg.fibers;
midpoints=zeros(length(allStreams),3);
endpointsRAS=midpoints;
endpointsLPI=endpointsRAS;
%speed it up a bit by using cell fun
streamSizes=cellfun(@size,allStreams,'UniformOutput',false);
streamLengths=cellfun(@(x) x(2),streamSizes,'UniformOutput',true);
midNodeNum=floor(streamLengths/2);

for iFibers=1:length(allStreams)
    curStreamline=allStreams{iFibers};
    midpoints(iFibers,:)=curStreamline(:,midNodeNum(iFibers));
    endpointsRAS(iFibers,:)=curStreamline(:,1);
    endpointsLPI(iFibers,:)=curStreamline(:,end);
end
% if a coordinate vectore is input, make sure it is rotated correctly.
% I.e. vertical column.

for iCriteriaInputs=1:length(criteriaSequence)-1

    currentCoordSpec=criteriaSequence{iCriteriaInputs};
    nextCoordSpec=criteriaSequence{iCriteriaInputs+1};
    
    switch lower(currentCoordSpec)
        case {'midpoint','midpoints'}
            currentCoord=midpoints;
        case 'endpointsras'
            currentCoord=endpointsRAS;
        case 'endpointslpi'
            currentCoord=endpointsRAS;
    end
    
   switch lower(nextCoordSpec)
        case {'midpoint','midpoints'}
            nextCoord=midpoints;
        case 'endpointsras'
            nextCoord=endpointsRAS;
        case 'endpointslpi'
            nextCoord=endpointsRAS;
    end
    
    
    %execute switch statement and find coords that meet criteria
    switch lower(relativePositioning)
        case 'superior'
            
            boolVec=currentCoord(:,3)>nextCoord(:,3);
            
        case 'inferior'
    
            boolVec=currentCoord(:,3)<nextCoord(:,3);
        case 'anterior'
     
            boolVec=currentCoord(:,2)>nextCoord(:,2);
        case 'posterior'
      
            boolVec=currentCoord(:,2)<nextCoord(:,2);
        case 'medial'
            if LRflag
                
                boolVec=currentCoord(:,1)>nextCoord(:,1);
            else
               
                boolVec=currentCoord(:,1)<nextCoord(:,1);
            end
        case 'lateral'
            if LRflag
               
                boolVec=currentCoord(:,1)<nextCoord(:,1);
            else
                
                boolVec=currentCoord(:,1)>nextCoord(:,1);
            end
            
    end
    
    % concat boolean vectors
    if iInputs==1
        booleanOut=boolVec;
    else
        booleanOut=and(booleanOut,boolVec);
    end
    
end
end