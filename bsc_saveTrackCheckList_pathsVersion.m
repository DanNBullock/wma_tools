function bsc_saveTrackCheckList_pathsVersion(csvPaths,plotProperties,zThresh,saveDir)


%workingDir='/N/dc2/projects/lifebid/HCP/Dan/EcogProject/proj-5c33a141836af601cc85858d'
%identifierTag='cleaned'
%plotProperties=[19 20]
%zThresh=[-4 6]

[subjNames,domainNames,propertyNames,valueArray]= bsc_normalizeStatMeasures_pathsVersion(csvPaths);

fileID = fopen(fullfile(saveDir,'checkList.txt'),'w');

fprintf(fileID,'\n Check performed on %s using thresh %s ',datestr(datetime),num2str(zThresh))

tractNames={domainNames{2:end}};
for iplotProperties=1:length(plotProperties)
    
    fprintf(fileID,'\n\n\n initiating check for %s',propertyNames{(plotProperties(iplotProperties)+1)})
    %indexing at 2 to get rid of wbfg
    plotArray=squeeze(valueArray(2:end,plotProperties(iplotProperties)+1,:));
    
    if length(zThresh)==2
        minz=min(zThresh);
        maxz=max(zThresh);
        [rowInd,colInd] = find(or((plotArray)>maxz,(plotArray)<minz));
    else
    [rowInd,colInd] = find(abs(plotArray)>zThresh);
    end
    [nanRowInd,nanColInd]=find(isnan(plotArray));
    
    
    troubleSubjects=unique([colInd' nanColInd']);

    for icheckSubjects=1:length(troubleSubjects)
       
        fprintf(fileID,'\n\n%s',subjNames{troubleSubjects(icheckSubjects)});
        checkvec=[];
        checkTheseTracts=rowInd(find(colInd==troubleSubjects(icheckSubjects)));
        if ~isempty(checkTheseTracts)
            for iTroubleTracts= 1:length(checkTheseTracts)
                checkvec=horzcat(checkvec,[tractNames{checkTheseTracts(iTroubleTracts)},' ',num2str(plotArray(checkTheseTracts(iTroubleTracts),troubleSubjects(icheckSubjects))),',']);
            end
            
            fprintf(fileID,'\n outliers: %s',checkvec);
        end
        clear checkTheseTracts
        
        emptyVec=[];
        emptyTracts=nanRowInd(find(nanColInd==troubleSubjects(icheckSubjects)));
        if ~isempty(emptyTracts)
            for iEmptyTracts= 1:length(emptyTracts)
                emptyVec=horzcat(emptyVec,[tractNames{emptyTracts(iEmptyTracts)},' ']);
            end
            if and(length(emptyTracts)==length(tractNames),~length(emptyTracts)==0)
                fprintf(fileID,'\n DATA MISSING');
            else
                fprintf(fileID,'\n failure: %s',emptyVec);
            end
        end
        clear emptyTracts    
    end
end
fclose(fileID)
end