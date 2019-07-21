function [avgTable, stdTable]=bsc_tableAverages_v3(csvPaths)

allProperties=[];
allDomains=[];

catDataStruc=[];
for isubjects =1:length(csvPaths)
    if exist(csvPaths{isubjects},'file')
        
        currTable=readtable(csvPaths{isubjects});
      
        
        currDomains=currTable{1:end,1};
        currProperties=currTable.Properties.VariableNames;
        %delete TractName
        namesIndex=find(~strcmp('TractName',currProperties));
        currData=currTable{1:end,namesIndex};
        currProperties={currProperties{namesIndex}};
        
        
        allProperties=unique(horzcat(allProperties,currProperties),'stable');
        %I guess they are vertical
        allDomains=unique(vertcat(allDomains,currDomains),'stable');
        
        newDataStruc=cell(length(allDomains),length(allProperties));
        
        curCatSize=size(catDataStruc);
        
        if ~isempty(catDataStruc)
            for iRows=1:curCatSize(1)
                for iCols=1:curCatSize(2)
                    newDataStruc{iRows,iCols}=catDataStruc{iRows,iCols};
                end
            end
        else
            %do nothing?
        end
        
        catDataStruc=newDataStruc;
        
         for iRows=1:length(currDomains)
                for iCols=1:length(currProperties)
                    rowInd=find(strcmp(currDomains{iRows},allDomains));
                    colInd=find(strcmp(currProperties{iCols},allProperties));
                    
                     %if not empty or isnan
                    if and(~isnan(currData(iRows,iCols)),~isempty(currData(iRows,iCols)))
                    newDataStruc{rowInd,colInd}=horzcat(newDataStruc{rowInd,colInd},currData(iRows,iCols));
                    else
                        %dont cat it
                    end
                    
                end
         end
    else
        %this case should never occur
        
    end
    catDataStruc=newDataStruc;
    %strucs should be done
end
        
avgTable=cellfun(@mean,catDataStruc);

stdTable=cellfun(@std,catDataStruc);
 
end