function [classificationOut] =bsc_streamlineCategoryPriors_TableBased(wbfg,atlas,inflateITer)
%  [classificationOut] =bsc_streamlineCategoryPriors_TableBased(wbfg,atlas,inflateITer)

%
% This function automatedly segments gross anatomical categories

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
%  classificationOut:  standardly constructed classification structure
%  Same for the other tracts
% (C) Daniel Bullock, 2019, Indiana University

[superficialClassification] =bsc_segmentSuperficialFibers(wbfg, atlas);
disp('superficial fibers identified')

%initialize classification structure
classificationOut=[];
classificationOut.names=[];
classificationOut.index=zeros(length(wbfg.fibers),1);

classificationMid=classificationOut;

%get the path to the repo
funcPath=which('bsc_streamlineCategoryPriors_TableBased');
repoPath=funcPath(1:strfind(funcPath,'wma_tools')-1);
grossAnatTablePath=fullfile(repoPath,'wma_tools','Utils','Data','GrossAnatomyLookup.csv');
grossAnatTable=readtable(grossAnatTablePath);

grossLabels=unique(grossAnatTable.GrossAnat);
hemiLabels=unique(grossAnatTable.Hemi);
origLabels=unique(grossAnatTable.x_No_);

endpoints1=zeros(3,length(wbfg.fibers));
endpoints2=zeros(3,length(wbfg.fibers));


for icategories=1:length(wbfg.fibers)
    curStream=wbfg.fibers{icategories};
    endpoints1(:,icategories)=curStream(:,1);
    endpoints2(:,icategories)=curStream(:,end);
end

disp('endpoints extracted')

if inflateITer>0
    [inflatedAtlas] =bsc_inflateLabels(atlas,inflateITer);
else
    inflatedAtlas=atlas;
end

atlasData=atlas.data;
grossAnatAtlas=atlas;
hemiAtlas=atlas;
grossAnatAtlasData=atlasData;
hemiAtlasData=atlasData;

%recode the atlases
for iOrigLabels=1:length(origLabels)
    currentLabel=origLabels(iOrigLabels);
    currentGrossLabel=find(strcmp(grossAnatTable{grossAnatTable.x_No_==currentLabel,{'GrossAnat'}},grossLabels));
    currentHemiLabel=find(strcmp(grossAnatTable{grossAnatTable.x_No_==currentLabel,{'Hemi'}},hemiLabels));
    grossAnatAtlasData(grossAnatAtlasData==currentLabel)=currentGrossLabel;
    hemiAtlasData(hemiAtlasData==currentLabel)=currentHemiLabel;
end
grossAnatAtlas.data=grossAnatAtlasData;
hemiAtlas.data=hemiAtlasData;

disp('atlases relabeled')

[endpoints1GrossIdentity] =bsc_atlasROINumsFromCoords_v3(grossAnatAtlas,endpoints1,'acpc');
[endpoints2GrossIdentity] =bsc_atlasROINumsFromCoords_v3(grossAnatAtlas,endpoints2,'acpc');
[endpoints1HemiIdentity] =bsc_atlasROINumsFromCoords_v3(hemiAtlas,endpoints1,'acpc');
[endpoints2HemiIdentity] =bsc_atlasROINumsFromCoords_v3(hemiAtlas,endpoints2,'acpc');

sameHemiBool=endpoints1HemiIdentity==endpoints2HemiIdentity;
%under the assumption that only groups 2 and 3 are left and right
%hemispheres
interHemiBool=[endpoints1HemiIdentity+endpoints2HemiIdentity]==5;
uFibersBool=superficialClassification.index>0;

disp('endpoint identities determined')

%START HERE TO PARSE CATEGORIES

%create blnk structure
identityStruc=cell(length(wbfg.fibers),1);
for iStreams=1:length(wbfg.fibers)
    %arbitrarily using 8 and 9 to indicate interhemi and ufiber,
    %respectiviely 
    %SPECIAL NOTE:  because we dont care about the ordering of the
    %endpoints, we can go ahead and sort them here, which will preclude us
    %having to deal with ordering issues later
    [~,~,curIdenties]=find([sort([endpoints1GrossIdentity(iStreams), endpoints2GrossIdentity(iStreams)]),interHemiBool(iStreams)*8,uFibersBool(iStreams)*9]);
    %assign identities to structure
    identityStruc{iStreams}=curIdenties;
end
%figure out what the unique names are.  Has to be done with strings because
%matlab?
streamNamesAsString=cellfun(@num2str,identityStruc,'UniformOutput',false);

%get the unique from these
uniqueIdentityCombinations=unique(streamNamesAsString);

%back convert into numbers so we can index into specific numeic elements
uniqueNamesAsNumbers=cellfun(@str2num,uniqueIdentityCombinations,'UniformOutput',false);

%can't figure out elegant way to do this.  I feel like there are decent
%ways to do this in numpy.  Regardless, here we are creating the identity
%and names vectors of the output classification.
outIdentityVec=zeros(length(wbfg.fibers),1);
namesVec=[];
for iIdentities=1:length(uniqueIdentityCombinations)
    curIdentityString=uniqueIdentityCombinations{iIdentities};
    curIdenetityNum=uniqueNamesAsNumbers{iIdentities};
    %vector containing the stream indexes associated with the current
    %identity
    currentStreamIdentityIndexes=find(ismember(streamNamesAsString,curIdentityString));
    outIdentityVec(currentStreamIdentityIndexes)=iIdentities;
    
    %generate current name
    currentName=strcat(grossLabels(curIdenetityNum(1)), '_to_',grossLabels(curIdenetityNum(2)));
    
    %append the special category if relevant.  
    %NOTE: in the current implementation ufiber and interhemi are mutually
    %exclusive, future versions of this would need more complex casing here
    %to extend this
    if length(curIdenetityNum)==3
        if curIdenetityNum(3)==8
            currentName=strcat(currentName,'_interHemi');
        elseif curIdenetityNum(3)==9
            currentName=strcat(currentName,'_uFiber');
        end
    end
    %enter current name into correct location
    namesVec{iIdentities}=currentName{1};
end

classificationOut.index=outIdentityVec;
classificationOut.names=namesVec;

classificationOut = wma_resortClassificationStruc(classificationOut);

disp('categorical segmentation complete')
end
