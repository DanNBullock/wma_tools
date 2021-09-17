for iImages=1:107
    convertFrame=num2str(iImages+1000);
    frameName=strcat('frame',convertFrame(2:end),'.png');
    sphereName=strcat('sphere_',num2str(iImages),'.png');
    
    niftiImage=imread(strcat('/N/u/dnbulloc/Carbonate/Downloads/proj-5a9ee80853bd38003c02afa7/sub-100206/figs/dwiFrames/',frameName));
    sphereImage=imread(strcat('/N/u/dnbulloc/Carbonate/Downloads/proj-5a9ee80853bd38003c02afa7/sub-100206/figs/spheres/',sphereName));
    
    %find Borders
    whiteMask=or(sum(sphereImage,3)==765,sum(sphereImage,3)==114);
    blackMask=sum(niftiImage,3)==0;
    
    [whiteCoordsX, whiteCoordsY]=find(~whiteMask);
    [blackCoordsX, blackCoordsY]=find(~blackMask);
    
    niftiXborderLeft=min(blackCoordsX);
    niftiXborderRight=max(blackCoordsX);
    
    niftiYborderBottom=min(blackCoordsY);
    niftiYborderTop=max(blackCoordsY);
    
    sphereXborderLeft=min(whiteCoordsX);
    sphereXborderRight=max(whiteCoordsX);
    
    sphereYborderBottom=min(whiteCoordsY);
    sphereYborderTop=max(whiteCoordsY);
    
    CroppedSphereData=sphereImage(sphereXborderLeft:sphereXborderRight,sphereYborderBottom:sphereYborderTop,:);
    shrunkSphere=imresize(CroppedSphereData,.35);
    
    sphereDims=size(shrunkSphere);
    
    %CroppedDWIData=niftiImage(niftiXborderLeft:niftiXborderRight+sphereDims(1),niftiYborderBottom:niftiYborderTop+sphereDims(2),:)
    CroppedDWIData=niftiImage(niftiXborderLeft:niftiXborderRight,niftiYborderBottom:niftiYborderTop,:);
    CroppedDWIData(end-sphereDims(1)+1:end,end-sphereDims(2)+1:end,:)=shrunkSphere;
    
    outFileName=strcat('/N/u/dnbulloc/Carbonate/Downloads/proj-5a9ee80853bd38003c02afa7/sub-100206/figs/croppedOverlay/fig_',num2str(iImages),'.png');
    imwrite(CroppedDWIData,outFileName)
end