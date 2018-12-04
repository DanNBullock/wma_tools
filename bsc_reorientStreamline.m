function [fiberTract] = bsc_reorientStreamline(fiberTract)
 %Inputs
 %fiberTract : a 3 x n vector, where n is the number of nodes in a
 %streamline.
 
 xDisplacement=abs(fiberTract(1,1)-fiberTract(1,end));
 yDisplacement=abs(fiberTract(2,1)-fiberTract(2,end));
 zDisplacement=abs(fiberTract(3,1)-fiberTract(3,end));
 
 displaceVec=[xDisplacement yDisplacement zDisplacement];
 maxDimIndex=find(displaceVec==max(displaceVec));
 
 if fiberTract(maxDimIndex,1) < fiberTract(maxDimIndex,end)
     fiberTract=fliplr(fiberTract);
 end
     
 