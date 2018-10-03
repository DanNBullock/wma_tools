%% this is a dummy script for collecting functions for the wma_seg_toolkit

mergedROI = wma_roiFromAtlasNum(pathToAtlas,ROInums,smoothKernel);

asegFile = wma_getAsegFile(fsDir , asegOption);

mergedROI = wma_roiFromFSnums(fsDir,fsROInums, smoothFlag, smoothKernel);

[roiOUT]=bsc_shiftROI(roiIN, shiftDir, shiftMag)

[fg, keep]=  bsc_tractByEndpointROIs(fg, rois)

bsc_plotEndpointsOnFG(fg)

[booleanOut]=bsc_applyMidpointCriteria(midpointsIn, varargin)

tractStruc = bsc_makeFGsFromClassification_v2(classification, wbFG,coordScheme)

[singleClassOUT]=bsc_extractClassification(classificationIN,nameOrIndex)

[classificationOut] = bsc_splitTractAtPoint(fg, coordinate,location ,classification)

[roiOUT] = bsc_planeFromROI(roiIN, location,fsDir)

[planarROI]=bsc_makePlanarROI_v2(fsDir,coord, dimension)

[classification] = bsc_reconcileClassifications(baseClassification,classificationAdd)

[classification] = bsc_spliceClassifications(baseClassification,classificationAdd)

[classificationOUT]=bsc_concatClassificationCriteria(classificationIN,name,varargin)

[roiOUT]=bsc_modifyROI(fsDir,roiIN, refCoord, location)

wma_SegmentFascicleFromConnectome(fg, rois, operation, fascicleFileName)

[classification] =wma_segLobes_v2(feORwbfg, fsDir)

[classification] =bsc_segmentAslant(wbfg, fsDir)

[classificationOut] =bsc_segmentMdLF_ILF_v3(wbfg, fsDir)

[classificationOut] = bsc_segpArcTPC(wbfg, fsDir)
