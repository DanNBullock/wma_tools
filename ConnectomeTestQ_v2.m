function [costFuncVec, AsymRat,FullDisp ,streamLengths, efficiencyRat ] = ConnectomeTestQ_v2(wbfg)

[wbfg, ~] = bsc_LoadAndParseFiberStructure(wbfg);

firstHalf=wbfg;
secondHalf=wbfg;
totalFib=length(wbfg.fibers);
totalFibPct=totalFib*.05;
for iFibers=1:totalFib
    if rem(iFibers,totalFibPct)==0
        fprintf('\n %i percent complete', (iFibers/totalFib)*100)
    end
    fiberNodeNum=round(length(wbfg.fibers{iFibers})/2);
    curStreamline=wbfg.fibers{iFibers};
    midpoints(iFibers,:)=curStreamline(:,fiberNodeNum);
    firstIndex=[1:fiberNodeNum];
    secondIndex=[fiberNodeNum:length(wbfg.fibers{iFibers})];
    firstHalf.fibers{iFibers}=curStreamline(:,firstIndex);
    secondHalf.fibers{iFibers}=curStreamline(:,secondIndex);
    endpoints1(iFibers,:)=curStreamline(:,1);
    endpoints2(iFibers,:)=curStreamline(:,end);
    firstHalfstreamLengths(iFibers)=sum(sqrt(sum(diff(firstHalf.fibers{iFibers},1,2).^2)));
    secondHalfstreamLengths(iFibers)=sum(sqrt(sum(diff(secondHalf.fibers{iFibers},1,2).^2)));
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
