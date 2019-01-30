function [costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ] = ConnectomeTestQ_v2(wbfg)

[wbfg, ~] = bsc_LoadAndParseFiberStructure(wbfg);


totalFib=length(wbfg.fibers);
totalFibPct=totalFib*.05;
allStreams=wbfg.fibers;
firstHalf=allStreams;
secondHalf=allStreams;
midpoints=zeros(length(allStreams),3);
endpoints1=midpoints;
endpoints2=midpoints;

firstHalfstreamLengths=midpoints(:,1);
secondHalfstreamLengths=firstHalfstreamLengths;
streamLengths=firstHalfstreamLengths;
firstDisp=firstHalfstreamLengths;
secondDisp=firstHalfstreamLengths;
firstDispRat=firstHalfstreamLengths;
secondDispRat=firstHalfstreamLengths;
FullDisp=firstHalfstreamLengths;
efficiencyRat=firstHalfstreamLengths;
AsymRat=firstHalfstreamLengths;

for iFibers=1:totalFib
    if rem(iFibers,totalFibPct)==0
        fprintf('\n %i percent complete', (iFibers/totalFib)*100)
    end
    fiberNodeNum=round(length(allStreams{iFibers})/2);
    curStreamline=allStreams{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    firstIndex=[1:fiberNodeNum];
    secondIndex=[fiberNodeNum:length(allStreams{iFibers})];
    firstHalf{iFibers}=curStreamline(:,firstIndex);
    secondHalf{iFibers}=curStreamline(:,secondIndex);
    endpoints1(iFibers,:)=curStreamline(:,1);
    endpoints2(iFibers,:)=curStreamline(:,end);
    firstHalfstreamLengths(iFibers)=sum(sqrt(sum(diff(firstHalf{iFibers},1,2).^2)));
    secondHalfstreamLengths(iFibers)=sum(sqrt(sum(diff(secondHalf{iFibers},1,2).^2)));
    streamLengths(iFibers)= firstHalfstreamLengths(iFibers)+secondHalfstreamLengths(iFibers);
    firstDisp(iFibers)=(sum(sum(squareform(pdist(vertcat(endpoints1(iFibers,:),midpoints(iFibers,:))))))/2);
    secondDisp(iFibers)=(sum(sum(squareform(pdist(vertcat(endpoints2(iFibers,:),midpoints(iFibers,:))))))/2);
    firstDispRat(iFibers)=firstDisp(iFibers)/firstHalfstreamLengths(iFibers);
    secondDispRat(iFibers)=secondDisp(iFibers)/secondHalfstreamLengths(iFibers);
    
    FullDisp(iFibers)=(sum(sum(squareform(pdist(vertcat(endpoints1(iFibers,:),endpoints2(iFibers,:))))))/2);
    efficiencyRat(iFibers)=FullDisp(iFibers)/streamLengths(iFibers);
    %imbalRat(iFibers)=(midEndRat1(iFibers)-midEndRat2(iFibers))/(midEndRat1(iFibers)+midEntRat2(iFibers));
    
    AsymRat(iFibers)=(firstDispRat(iFibers)-secondDispRat(iFibers))^2;
    
    
    %midEndRat1(iFibers)=(sum(sum(squareform(pdist(vertcat(endpoints1(iFibers,:),midpoints(iFibers,:))))))/2)/streamLengths(iFibers);
    %midEndRat2(iFibers)=(sum(sum(squareform(pdist(vertcat(endpoints2(iFibers,:),midpoints(iFibers,:))))))/2)/streamLengths(iFibers);
end




%costFuncVec=times((1-AsymRat).^-1,streamLengths);
costFuncVec=(1-AsymRat).^-1;
end
