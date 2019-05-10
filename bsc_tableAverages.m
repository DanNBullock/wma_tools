function [avgTable, stdTable]=bsc_tableAverages(csvPaths)
catDomains=[];
catData=[];
for isubjects =1:length(csvPaths)
    if exist(csvPaths{isubjects},'file')
        currTable=readtable(csvPaths{isubjects});
        
        currDomains=currTable{1:end,1};
        currData=currTable{1:end,2:end};
        
        %if your first subject is empty given your criteria its going to
        %cause a problem.  Also you need to sort out your life.
        if isempty(catData)
            catData=currData;
            propertyNames=currTable.Properties.VariableNames;
            avgTable=currTable;
            stdTable=currTable;
            
            
        else
            
            %checks to see that same domainsize is being merged
            checkCat=length(catDomains)==length(currDomains);
            if checkCat
                %if they are the same size, we still need to check to see
                %if the order of the items matches
                [~,ia,ib] = intersect(catDomains,currDomains,'stable');
                if~isequal(ia,ib)
                    %could error here if you have two diffrent
                    %classifications with exactly the same number of items
                    %in them
                     addData=[];
                    for iDomains=1:length(currDomains)
                        curIndx=find(strcmp(currDomains{iDomains},catDomains));
                        if ~isempty(curIndx)
                        spliceData(curIndx,:)=currData(iDomains,:);
                        else
                        addData=vertcat(addData,currData(iDomains,:));
                        end
                    end
                     spliceData=vertcat(spliceData,addData);
                    %once you've rebuilt the structure transfer it to
                    %currData
                    currData=spliceData;
                    clear spliceData
                else
                    %if they are equal, then nothing needs to be done
                end
                %if they are of different lengths then we now have to do
                %some work to figure out what's going on
            else
                if length (catDomains)>length(currDomains)
                [diffDomains] = setdiff(catDomains,currDomains,'stable');
                [~, ia, ib]=intersect(catDomains,diffDomains,'stable');
                [~,presentA,presentB] = intersect(catDomains,currDomains,'stable');
                    for iDomains=1:length(catDomains)
                        if ~ismember(iDomains,ia)
                            %might screw up if iDomains isn't tantamount to
                            %the presentA, but whatever?
                        correspondingInd=find(iDomains==presentA); 
                        spliceData(iDomains,:)=currData(correspondingInd,:);
                        else
                        spliceData(iDomains,:)=NaN(1,length(catData(1,:,1)));
                        end
                    end
                    currData=spliceData;
                    clear spliceData
                elseif length (catDomains)<length(currDomains)
                    addData=[];
                    %maybe necessary, but not sued currently.  Possibly the
                    %way used is inefficeint
                    %[diffDomains] = setdiff(currDomains,catDomains,'stable');
                    %[~, ia, ib]=intersect(currDomains,diffDomains,'stable');
                    %[~,presentA,presentB] = intersect(currDomains,catDomains,'stable');
                    for iDomains=1:length(currDomains)
                        %probably a better way to do this
                        curIndx=find(strcmp(currDomains{iDomains},catDomains));
                        if ~isempty(curIndx)
                        spliceData(curIndx,:)=currData(iDomains,:);
                        else
                        addData=vertcat(addData,currData(iDomains,:));
                        end
                    end
                    currData=vertcat(spliceData,addData);
                %deal with it when it comes up #lazy
                
                end 
            end
            %now we cat the data, however we must first make sure they are
            %the right dimesnions
            
            %this could trigger if the catData were larger than the curr
            %data, but the earlier sections should have dealt with that
            catSize=size(catData);
            curSize=size(currData);
            if ~isequal(catSize,curSize)
                %if you catData isn't actually 3d yet, compensate for that
                if length(catSize)==3
                    nanAppend=NaN(curSize(1)-catSize(1),catSize(2),catSize(3));
                else
                    nanAppend=NaN(curSize(1)-catSize(1),catSize(2));
                end
                
                catData=cat(1,catData,nanAppend);
            else
           %do nothing
            end
              catData=cat(3,catData,currData);
        end
    catDomains=cat(1,currDomains,catDomains);
    catDomains=unique(catDomains,'stable');
    
    else
        currentCatDim=size(catData);
        %lol
        nanCat=nan(currentCatDim(1),currentCatDim(2));
        catData=cat(3,catData,nanCat); 
    end    
end

stdData=std(catData,0 ,3,'omitnan');
meanData=mean(catData,3,'omitnan');

dataSize=size(meanData);
%varTypes=cell(1,dataSize(2)+1);
varTypes{1}='string';
for iProperties=1:dataSize(2)
varTypes{1+iProperties}='double';
end
%avgTable= table('Size',[dataSize(1),dataSize(2)+1],'VariableTypes',varTypes);
%stdTable= table('Size',[dataSize(1),dataSize(2)+1],'VariableTypes',varTypes);

[DataRows,DataColums]=size(meanData);

avgHold=cell(DataRows,DataColums+1);
avgHold{1:end,1}=catDomains;
avgHold{1:end,2:end}=meanData;

avgTable=cell2table(avgHold);
avgTable.Properties.VariableNames=propertyNames;

stdHold=cell(DataRows,DataColums+1);
stdHold{1:end,1}=catDomains;
stdHold{1:end,2:end}=stdData;

stdTable=cell2table(stdHold);
stdTable.Properties.VariableNames=propertyNames;

end