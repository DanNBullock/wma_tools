function [allDomains,allProperties,valueArray]= bsc_normalizeStatMeasures_pathsVersion_v2(csvPaths)

%workingDir='/N/dc2/projects/lifebid/HCP/Dan/EcogProject/proj-5c33a141836af601cc85858d'
%identifierTag='measures'

%csvPaths = tractStatNamesGen(workingDir,identifierTag);

[avgTable, stdTable]=bsc_tableAverages_v3(csvPaths);


%stolen almost verbatum from bsc_tableAverages_v3
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
         
         curLengths=cellfun(@length,newDataStruc);
         maxLength=max(max(curLengths));
         [nonMaxRowInds,nonMaxColInds]=find(~[curLengths==maxLength]);
         for iNonMax=1:length(nonMaxRowInds)
         newDataStruc{nonMaxRowInds(iNonMax),nonMaxColInds(iNonMax)}=horzcat(newDataStruc{nonMaxRowInds(iNonMax),nonMaxColInds(iNonMax)},nan);
         end
         clear nonMaxRowInds
         clear nonMaxColInds
         
    else
        %this case should never occur
        
    end
    catDataStruc=newDataStruc;
    %strucs should be done
end

%there's a smarter way to do this.  At least it doesn't take long
finalCatSize=size(catDataStruc);

valueArray=[];
 for iRows=1:finalCatSize(1)
                for iCols=1:finalCatSize(2)
                 valueArray(iRows,iCols,:)= catDataStruc{iRows,iCols};
                end
 end

finalDimensions=size(valueArray);
 
for isubjects=1:finalDimensions(3)
    valueArray(:,:,isubjects)=rdivide(minus(valueArray(:,:,isubjects),avgTable),stdTable);  
end
     
end