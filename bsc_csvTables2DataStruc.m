function [allProperties, allDomains, outDataStruc]=bsc_csvTables2DataStruc(csvPaths)
%function [allProperties, allDomains, dataStruc ]=bsc_csvTables2DataStruc(csvPaths)
%
%  Purpose:  This function reads in a list of csv tables (ostensibly
%  covering the same type of table) and outputs a single data structure. 
%  Designed to take care of both situations where not some tables may have
%  more or fewer domains (things that can have properties).  Probably wont
%  handle more or fewer properties.  That wouldn't happen under normal
%  usage cases.  Also, this whole thing probably only works for
%  numerically based tables.
%
%  INPUTS
%  csvPaths:  a cell strucure containing some number of paths to csv tables
% 
%  OUTPUTS
%
%  allProperties:  A list of all properties covered in the output data
%  structure.  Think of it as the column labels.  
%
%  allDomains:  A list of all domains covered in the output data
%  structure.  Think of it as the row labels.  
%
%  dataStruc:  The output data structure itself
%
%  Created by Dan Bullock 5-27-2019
%%  Begin code

%initialize data structures
allProperties=[];
allDomains=[];
catDataStruc=[];

for isubjects =1:length(csvPaths)
    %try to read the file if it exists
    if exist(csvPaths{isubjects},'file')
        
        %read the table
        currTable=readtable(csvPaths{isubjects});
        %get the property/variable names from the table properties
        currProperties=currTable.Properties.VariableNames;
        %check to see if domains are being specified by the table labels or
        %by the first column.
        if isempty(currTable.Properties.RowNames)
            currDomains=currTable{1:end,1};
            %lets assume that the first column is somehow specifying the
            %domains, as such we should eliminate the first column from
            %subsequent steps
            namesIndex=2:length(currProperties);
            currData=currTable{1:end,namesIndex};
            currProperties={currProperties{namesIndex}};
        else
            currDomains=RowNames;
        end

        % Horcat the property names together.  Given that we are assuming
        % that properties are necessarily shared across tables (i.e. no
        % table type mixing) this is probably unnecessary.
        allProperties=unique(horzcat(allProperties,currProperties),'stable');
        
        %I guess they are vertical, it is the row labels, after all.
        allDomains=unique(vertcat(allDomains,currDomains),'stable');
        
        %now comes the trick.  Create a blank cell structure THAT IS THE
        %SIZE OF THE PREVIOUS ITERATION'S OUTPUT, and then fill it with nan
        %vecs that are the isubject-1 in length.  
        newDataStruc=cell(length(allDomains),length(allProperties));
        nanVec=NaN(1,length(isubjects)-1);
        %apply nan to entire struc
        newDataSize=size(newDataStruc);
        for iRows=1:newDataSize(1)
            for iCols=1:newDataSize(2)
                newDataStruc{iRows,iCols}=nanVec;
            end
        end
        
        %get the size of the catDataStruct
        curCatSize=size(catDataStruc);
        %if it is not empty, replace the nanVecs with the values obtained
        %from all previous tables.  In the event that the current table has
        %more domains than previous ones, this will leave the nan vec (of
        %isubjects-1 length) in place so that this subjects data can be
        %concatonated without problems, and the nan values can be ascribed
        %to previous subjects which had no values for that domain.
        if ~isempty(catDataStruc)
            for iRows=1:curCatSize(1)
                for iCols=1:curCatSize(2)
                    newDataStruc{iRows,iCols}=catDataStruc{iRows,iCols};
                end
            end
        else
            %do nothing?
        end
        
        %now we append the data specific to this subject
         for iRows=1:length(currDomains)
                for iCols=1:length(currProperties)
                    %we assume nothing about the ordering of the properties
                    %or the domains, instead, we find the correspondnace of
                    %the current table's index to the newDataStruc's (the
                    %amalgamated structure we are building) index.  This
                    %presumably improves robustness.
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
         
         %in the event that a value is not appended to the relevant cell of
         %newDataStruc.  Its length will be one less than the current
         %subject number.  Or more accurately, the current max vector
         %length.  It could be the case that a bad path to a csv path is
         %input, in which case this code would skip that csv, and there
         %would now be a mismatch between the max vec length and the
         %subject number.
         
         % get the lengths of all cells data vectors
         curLengths=cellfun(@length,newDataStruc);
         %find the greatest amongst them (under most circumstances, this
         %should be ijsubjects)
         maxLength=max(max(curLengths));
         %find those indexes which do not have the appropriate number of
         %items in them
         [nonMaxRowInds,nonMaxColInds]=find(~[curLengths==maxLength]);
         %append nans when apprpopraite.  Takes care of nans, empty values
         %and absent domains.
         for iNonMax=1:length(nonMaxRowInds)
         newDataStruc{nonMaxRowInds(iNonMax),nonMaxColInds(iNonMax)}=horzcat(newDataStruc{nonMaxRowInds(iNonMax),nonMaxColInds(iNonMax)},nan);
         end
         %clear the vectors for the next iteration.  Probably not
         %necessary.
         clear nonMaxRowInds
         clear nonMaxColInds
    else
        %this case should never occur, but throw a warning anyways
        warning('csv %s not found, skipping',csvPaths{isubjects})
    end
    %shift the newDataStruc to the catDataStruc holder
    catDataStruc=newDataStruc; 
end
% catDataStruc should be complete and represent entire data set now.
finalCatSize=size(catDataStruc);
%covert current cell structure to a standard 3d array.  May not have been
%necessary to do it this way, but oh well.
outDataStruc=zeros(finalCatSize(1), finalCatSize(2),length(catDataStruc{1,1}));
for iRows=1:finalCatSize(1)
    for iCols=1:finalCatSize(2)
        outDataStruc(iRows,iCols,:)= catDataStruc{iRows,iCols};
    end
end
%done
end