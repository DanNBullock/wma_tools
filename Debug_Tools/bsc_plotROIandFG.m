function bsc_plotROIandFG(fg,roi,Incolor)
%function bsc_plotROIandFG(fg,roi,Incolor)
%

if ischar(roi)
    roi=bsc_loadAndParseROI(roi);
end
   
bsc_plotEndpointsOnFG_v2(fg)

scatter3(roi.coords(:,1),roi.coords(:,2),roi.coords(:,3),3,Incolor)

end