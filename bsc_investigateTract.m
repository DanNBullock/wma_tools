function [endZone1, endZone2 ]=bsc_investigateTract(wbfg,atlas)

allStreams=wbfg.fibers;

%atlasPath=fullfile(fsDir,'/mri/','aparc.a2009s+aseg.nii.gz');

for icategories=1:length(allStreams)
    curStream=allStreams{icategories};
    endpoints1(:,icategories)=curStream(:,1);
    endpoints2(:,icategories)=curStream(:,end);
end

fprintf('endpoints extracted\n')

inflateITer=2;

if inflateITer>0
    [inflatedAtlas] =bsc_inflateLabels(atlas,inflateITer);
else
    inflatedAtlas=atlas;
end

[endZone1] =bsc_atlasROINumsFromCoords_v3(inflatedAtlas,endpoints1,'acpc');
[endZone2] =bsc_atlasROINumsFromCoords_v3(inflatedAtlas,endpoints2,'acpc');
fprintf('endpoint identities determined\n')

end
