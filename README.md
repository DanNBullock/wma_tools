

# White Matter Anatomy tools (wma_tools)

## Overview description

This repository contains a number of matlab code tools which can be used to perform anatomically informed segmentations and analyze the outputs.  Additionally, it includes several examples of _actually_ implemented segmentations for established white matter tracts.

#### Pythonic alternative availability
**March 9, 2020 notice**: Comperable pythonic implementations of most **wma_tools** functionalities can be found in the [wma_pyTools
](https://github.com/DanNBullock/wma_pyTools/tree/packageVers) repository.  The [wma_pyTools
](https://github.com/DanNBullock/wma_pyTools/tree/packageVers) repository also includes a host of [visualization](https://github.com/DanNBullock/wma_pyTools/blob/packageVers/wmaPyTools/visTools.py) and [analysis](https://github.com/DanNBullock/wma_pyTools/blob/packageVers/wmaPyTools/analysisTools.py) capabilties not inculded in the Matlab-based **wma_tools** toolkit.

## Author, funding sources, references

### Authors
- Daniel Bullock (dnbulloc@iu.edu)

### Laboratory Director
- Franco Pestilli (franpest@indiana.edu, pestilli@utexas.edu)

### Funding

The development of this code was directly supported by the following funding sources.

[![NIMH-T32-5T32MH103213-05](https://img.shields.io/badge/NIMH_T32-5T32MH103213--05-blue.svg)](https://projectreporter.nih.gov/project_info_description.cfm?aid=9725739)  

[![NICT-ITP-2016&17](https://img.shields.io/badge/NICT_ITP-2016&17-ITP-9cf.svg)](https://www.nict.go.jp/en/global/internship.html)

### References 

[Bullock, D., Takemura, H., Caiafa, C. F., Kitchell, L., McPherson, B., Caron, B., & Pestilli, F. (2019). Associative white matter connecting the dorsal and ventral posterior human cortex. Brain Structure and Function, 224(8), 2631-2660.](https://doi.org/10.1007/s00429-019-01907-8)

[Wassermann, D., Makris, N., Rathi, Y., Shenton, M., Kikinis, R., Kubicki, M., & Westin, C. F. (2016). The white matter query language: a novel approach for describing human white matter anatomy. Brain Structure and Function, 221(9), 4705-4721.](https://doi.org/10.1007/s00429-015-1179-4)

## Repository Overview

### Directory structure and brief descriptions

The following table provides an overview of the types of code tools contained within each directory of the repository

| Directory                 | Description                                                                                                                                       |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| [Analysis](#analysis)                  | Analysis tools for segmented tracts; both individual and group level                                                                              |
| [Atlas_tools](#atlas_tools)               | Tools for working with NIfTI atlases/parcellations; both processing and analysis                                                                  |
| [BL_Wrappers](#bl_wrappers)               | Wrappers for interacting with brainlife.io apps                                                                                                   |
| [ClassificationStruc_Tools](#classificationstruc_tools) | Tools for working with [white matter classification (WMC) structures](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b)                     |
| [Debug_Tools](#debug_tools)               | Tools for troubleshooting and developing segmentations; typically lightweight streamline visualizations.                                          |
| [ROI_Tools](#roi_tools)                 | Tools for obtaining, modifying, and utilizing ROIs.  Both NiFTI and [vistasoft](https://github.com/vistalab/vistasoft), point-cloud format.       |
| [Segmentations](#segmentations)             | Implemented segmentations using wma_tools methods; contains single tract and multi-tract implementations; contains archived segmentation versions |
| [Stream_Tools](#stream_Tools)              | Tools for assessing streamline characteristics, superficial modification, and criteria application                                                |
| [Utils](#utils)                     | General use utility functions, externally sourced functions                                                                                       |
| [Visualization](#visualization)             | Functions for generating visualizations (typically publication quality)                                                                           |

## Detailed (but summarized/non-exhaustive) directory content descriptions

Below, the sections will discuss the specific tools/functions contained within each directory.  Relevant data domain(s) will be designated below general-form function name.  Toolkit designations (i.e. "wma_[function name]", for white matter anatomy, and "bsc_[function name]", for bloomington script compilation) are dropped when designating function names.  "Version descrepancies/differences will not be discussed.  Higher numbered versions should be presumed to be the latest iteration of a tool/function.

### Analysis

This directory contains several functions/tools which are generally used to compute derived statistics or summaries from other data objects (i.e. CSVs and tractomes/[white matter classification (WMC) structures](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b))

#### ConnectomeTestQ
Relevant data domain(s):  tractography
A connectome quality test.  
A function which computes moderately obscure quantative features of streamlines from an input tract (or tractome).

| Output Characteristic/Var | Description                                                                                                                                       |
|---------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------|
| [streamLengths](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L42)             | The length traversed by each streamline                                                                              |
| [FullDisp](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L48)                  | The actual spatial [displacement](https://en.wikipedia.org/wiki/Displacement_(geometry)) of each streamline                                                                  |
| [efficiencyRat](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L49)             | Efficiency ratio: the ratio between the displacement and the actual traversal length of each streamline.  1 = maximal efficiency, 0 = minimal efficiency (e.g. ending up exactly where it started).  A perfectly semicircular streamline (i.e. u shaped) would have an efficiency ratio of 1/π.                                                                                                 |
| [AsymRat](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L52)                   | The square of the difference between the the **efficiencyRat** for the two halves of the streamline.  1 = maximal asymmetry, 0 = minimal asymmetry. <br/>  _Example 1_ : if the first half was completely circuitous and the second half was a straight line;  (0-1)^2 =1 <br/> _Example 2_ : both halves are a straight line; (1-1)^2 =0 ;                  |
| [costFuncVec](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L91)               | The inverse of 1-**AsymRat** .  As such, inf in cases of perfect asymmetry, 0 in cases of minmal asymetry (i.e. straight line).  Considered as a biological prior/cost function for efficient/sensible wiring.                                           |

#### csvTables2DataStruc (csv, general data)
Relevant data domain(s):  csv, general data

A function for merging multiple 2-D data tables into a 3-D data structure.  It is assumed that the individual csv tables are storing information in rougly identical layouts (i.e. subject data).  Robust against row ("Domain") and column ("Property") droppage (i.e. absence in some files).  Essentially designed to amalgamate subject level data (ostensibly stored in distinct csv files) into a group level data structure.  Possibly a hacky variant of a [pandas](https://pandas.pydata.org/) function/method.

#### feAndAFQqualityCheck
Relevant data domain(s):  tractography
 
A function for visualizing (for a single tractome/[LiFE input](https://brainlife.io/datatype/58d15eaee13a50849b258844)) **many** quantative tractography traits.  Adaptively alters analyses, quantative output, and figure layout depending on input.  If optional [LiFE structure](https://brainlife.io/datatype/58d15eaee13a50849b258844) and/or [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) are input, will generate statistics specific to those structures as well (i.e. "surviving streamlines" and "white matter tracts", respectively) and the combination thereof (surviving streamlines in segmented white matter tracts).  [quantAllWMNorm](#quantallwmnorm) and [quantWBFG](quantwbfg) constitute workhorse functions.  
See [this readme](https://github.com/brainlife/app-tractographyQualityCheck/tree/1.2) for an explicit listing, description, and code linkage of each quantative feature computed.  
Featured on Brainlife.io as a standalone app:  
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.189-blue.svg)](https://doi.org/10.25663/brainlife.app.189)

#### normalizeStatMeasures_pathsVersion
Relevant data domain(s):  csv, general data

Computes group-level normalized statistics for all rows ("Domains") and columns ("Properties") in a multi-subject (i.e. csv) dataset.  Cross-subject average computed with [tableAverages](#tableaverages).  Associated with functionality of [csvTables2DataStruc](csvtables2datattruc).  Likely a hacky version of a [pandas](https://pandas.pydata.org/) or SQL function/method.

#### plotZscoreMeasuresFromCSV
Relevant data domain(s):  csv, general data

Essentially (though not actually), computes, saves, and plots the Z-score normalized output of [normalizeStatMeasures_pathsVersion](#normalizestatmeasures_pathsversion) . Good for publication figure generation.  Optional inputs **plotProperties** and **subSelect** allow for selection of particular rows ("Domains") and columns ("Properties") (and combinations thereof) for plotting/computation.

#### streamlineGeometryPriors
Relevant data domain(s):  (tractography)

Runs the analysis performed by [ConnectomeTestQ](#connectometestq) in order to create a [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) for both the [_Asymmetry Ratio_](#connectometestq) and [_Efficiency Ratio_](#connectometestq) , wherein the distinct streamline identification categories contained within the [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) constitute the .05 incrimented bins of those quantative features (thus resulting in 20 distinct streamline categories per structure).  Intended to be used as a biological/geometric prior for segmentation algorithms.

#### bsc_tableAverages
Relevant data domain(s):  csv, general data

Workhorse function for computing averages across tables with roughly identical data layouts (i.e. subject data).  Robust against row ("Domain") and column ("Property") droppage (i.e. absence in some files).  Likely a hacky version of a [pandas](https://pandas.pydata.org/) or SQL function/method.

#### endpointMapsDecay
Relevant data domain(s):  tractography

Generates [smoothed/inflated](https://www.mathworks.com/help/matlab/ref/smooth3.html) NIfTI **density** (i.e. count) masks **for the streamline endpoints** of tracts identified in the input [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  Requires a reference NIfTI (to establish mask NIfTI dimensions) and a refrence tract/tractome (to extract identified streamlines from).  Relies on **endpointClusterProto** in order to ensure streamline endpoints are associated with appropriate endpoint group/NIfTI output.  Workhorse function for **[classifiedStreamEndpointCortex]**

#### quantAllWMNorm
Relevant data domain(s):  tractography

Quantifies more common quantative traits associated with white matter.  Underwrites [feAndAFQqualityCheck](#feandafqqualitycheck).  [quantAllWMNorm](#quantallwmnorm) is predicated upon the assumption of an associated [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) and computes a number of stats relative to the (presumed) source whole-brain tractome.  This is as opposed to [quantWBFG](#quantwbfg) which makes no presumptions regarding associated [white matter classification (WMC) structures](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  
See [this readme](https://github.com/brainlife/app-tractographyQualityCheck/tree/1.2) for an explicit listing, description, and code linkage of each quantative feature computed.

#### quantWBFG
Relevant data domain(s):  tractography

The standalone quantification function for a singleton (presumed whole-brain tractome) tract structure.  Quantifies more common quantative traits associated with white matter.  Associated with [quantAllWMNorm](#quantallwmnorm).  
See [this readme](https://github.com/brainlife/app-tractographyQualityCheck/tree/1.2) for an explicit listing, description, and code linkage of each quantative feature computed.

### Atlas_tools

This directory contains a number of tools and functions for working with NIfTI atlases/parcellations.

#### computeAtlasStats

This function iterates across the unique label values (presumed integer, as would be the case for an atlas; bad things would happen for continuous float NIfTI) found in a NIfTI and computes the following statistics

| Quantity label       | Description |
|----------------------|--------------------------------------------------------------------------------------------------------|
| actualVol            | The actual, world/subjectspace volume of the ROI        |
| wholeBrainProportion | The proportion of the total brain  (non 0 label entries) volume occupied by the ROI                          |
| centroidx            | x coordinate of centroid        |
| centroidy            | y coordinate of centroid      |
| centroidz            | z coordinate of centroid     |
| medialBorderCoord    | x coordinate of medial border       |
| lateralBorderCoord   | x coordinate of lateral border     |
| anteriorBorderCoord  | y coordinate of anterior border      |
| posteriorBorderCoord | y coordinate of posterior border      |
| superiorBorderCoord  | z coordinate of superior border       |
| inferiorBorderCoord  | z coordinate of inferior border     |
| boxyness             | the ratio of the actual ROI volume to the rectangular prism/cuboid formed by its borders.  A perfectly rectangular prism/cuboid ROI = 1, a perfectly spherical ROI = π/6 |

This function is featured as a component of an app on [brainlife.io](https://brainlife.io/):  
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.272-blue.svg)](https://doi.org/10.25663/brainlife.app.272)  
(although this app is specific to [freesurfer output](https://brainlife.io/datatype/58cb22c8e13a50849b25882e) this function can be run on any [atlas/parcellation](https://brainlife.io/datatype/5c1a7489f9109beac4a88a1f)

#### computeFresurferStats

A precursor to [computeAtlasStats](#computeatlasstats), which takes a [freesurfer directory](https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/OutputData_freeview) as input (and thereby only works on "aparc.a2009s+aseg.nii.gz", assuming it exists, which it doesn't by default, due to the default .mgz output) and computes a subset of the same quantative features (e.g. no actual volume, centroid, or boxyness).

#### inflateLabels

A freesurfer-specific, "aparc.a2009s+aseg.nii.gz" [grey](https://github.com/DanNBullock/wma_tools/blob/3bab47ab121c273e43e1f2fb48647c41e66bfb00/Atlas_tools/bsc_inflateLabels.m#L22)/[subcortical](https://github.com/DanNBullock/wma_tools/blob/3bab47ab121c273e43e1f2fb48647c41e66bfb00/Atlas_tools/bsc_inflateLabels.m#L23)-specific inflation function which inflates gray and subcortical labels into the white matter.  Useful for preprocessing the fs_a2009_aseg for subsequent use in segmentation, particularly in the case of tractography generation algorithms which do not guarentee termination in grey / subcortical structures (e.g. mrtrix2).  Relies on modified versions of [fnDeislandLabels](fndeislandlabels), originally written by [Brent McPherson](https://github.com/bcmcpher)
Could theoretically be adapted to work with other pacellations assuming the provision of distinct grey and white matter mask inputs.

#### inflateRelabelIslands

A "user ready" (i.e. no user speficication of parameters, with an important caviot), atlas-general (i.e. not specific to freesurfer) function which detects and relabels parcellation islands (label subsections which are not connected to the largest, contiguous label component).  Important caviot is that the input parcellation cannot feature the label 999, which is used as a [trash label](https://github.com/DanNBullock/wma_tools/blob/3bab47ab121c273e43e1f2fb48647c41e66bfb00/Atlas_tools/bsc_inflateRelabelIslands.m#L32) in this function.  [Could and should be fixed in future versions (i.e. -1)](https://github.com/DanNBullock/wma_tools/issues/21).

#### fnDeislandLabels

A modified version of the fnDeislandLabels code originally written by [Brent McPherson](https://github.com/bcmcpher), adapted for current author's norms.  essentially the same as [inflateRelabelIslands](#inflaterelabelislands), except that it doesn't specify the [maxisleSize](https://github.com/DanNBullock/wma_tools/blob/3bab47ab121c273e43e1f2fb48647c41e66bfb00/Atlas_tools/fnDeislandLabels_v3.m#L17) and ["trash label" (replaceVal)](https://github.com/DanNBullock/wma_tools/blob/3bab47ab121c273e43e1f2fb48647c41e66bfb00/Atlas_tools/fnDeislandLabels_v3.m#L20), which must be provided by the user.

#### getAsegFile

A utility function which loads one of the .aseg parcellation options ('orig' or '2009') output from [freesurfer](https://surfer.nmr.mgh.harvard.edu/).  If the .nii.gz version of the parcellation isn't found, this function will attempt to use [mri_convert](https://surfer.nmr.mgh.harvard.edu/fswiki/mri_convert).  As such, some of the adaptive functionality of this function is predicated upon freesurfer/freeview being installed on the local file structure.

### BL_Wrappers

The contents of the **BL_Wrappers** will not be discussed in this overview.  Previously, this directory contained wrappers that were needed to interface with earlier versions of the [brainlife.io] platform and/or generate appropraite output.  These functionalities are now typically handled by distinct scripts which depend on the datatypes involved.

### ClassificationStruc_Tools

The **ClassificationStruc_Tools** directory contains functions which are used for generating, modifying and/or interacting with [white matter classification (WMC) structures](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).

#### concatClassificationCriteria

This function essentially applies a logical _and_ operation across a set of input boolean vectors (each with length equal to the number of streamlines in the structure associated with the source tractogram) and creates a new entry in the input [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) using the specified name.  The presumption/use case for this is that the input boolean vectors represent a number of criteria that have been applied to a sourcetractogram such that the TRUE values are associated with streamlines that meet the criteria and FALSE values do not.  In this way, the resultant tract classification added to the classification structure represents those streamlines which meet all criteria.  This function exists because not all criteria used in segmentation are based on the standard ROI intersection method.

#### extractClassification

This function creates a new classification structure with the **only** tract labeled being the specified tract or index.  The use case for this is to avoid  having to pass around segmented tracts in order to modify them.  Instead you pass around the whole tractogram and the classification structure.  Allows for lightweight tract representation and non-duplication of streamlines

#### extractStreamIndByName

This function to a extract boolean vector indexing the the streamlines associated with the associated string label.  Used to improve robustness and facilitate error handling.  Outputs only a boolean vector, not a [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b)

#### makeFGsFromClassification

Extremely important/useful function which creates a stucture containing *all* of the identified/labeled streamline groups contained within a classification structure.  This facilitates secondary uses like postprocessing, plotting, or other visualization/analysis.  
NOTE: the output tractStruc is a structure wherein each row corresponds to a diffrent classification group. Previous versions of this function added a layer of complexity to this, such there were .fg and .name fields.  This was redundant because there is a .name field within .fg by default.  As such, this is now just a basic [cell structure](https://www.mathworks.com/help/matlab/cell-arrays.html) (wherein each cell contains a distinct "fg object"), and thus there are no fields.

#### makeWBFGandClassFromList

Conceptually, this function is the inverse of [makeFGsFromClassification](makefgsfromclassification).  Instead of generating separate tract structures/files from an input tractome and [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b), this function takes in a list of [.tck format files](https://brainlife.io/datatype/5dcf0047c4ae28d7f2298f48) and combines them into a single tract structure ([vistasoft](https://github.com/vistalab/vistasoft) fg format) and also creates a corresponding classification structure which identifies the inherited tract identies of the new structure's streamlines.  Used to convert from use cases wherein distinct tracts are stored as separate file objects.

#### mergeFGandClass

This function merges input fiber goups and classification structures.  A deceptively complicated function which relies on [reconcileClassifications](#reconcileclassifications), [spliceClassifications](#spliceclassifications) to handle subtly distinct use cases.   
Example use cases:

- 1.  Merging a tracking parameter sweep:  If you have several different fiber groups generated by different tractography parameter settings, this will provide a merged fg structure along with a classification structure indicating which streamlines came from which source tractography method (as specified by the file name input)

-  2.  Merging multiple segmentations:  By inputting multiple source whole brain fiber groups (from the same subject, one assumes), along with the classification structure from each corresponding segmentation, you can create a "merged replication segmentation" which presumably has more streamlines representing each
 particular tract, or better represents distinct sets of tracts (e.g. if tailored source tractography methods were used).
 
-  3.  Incorporating tracts generated from other methods into an exsting whole brain fiber group and segmentation:  if one already has a whole brain tractome (and corresponding segmentation) it is possible to "append" more tracts to the source tractome and classification structure.  The tracts to be appended *do not need corresponding classification structures* and will be added to the existing classification structure (the one that corresponds to the whole brain tractography, should it exist).  Their names will be derived from the name in the fg.name field.

#### reconcileClassifications

This function is for reconciling two [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  In reconciling, the presumption is that these correspond to the **same** source tractography entities (i.e. collections of streamlines / source tractogram) _but the entities classified be the WMC structure could be the same or different in the two input structures_.

#### spliceClassifications

This function is for splicing two [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  When splicing, the presumption is that these correspond to **separate** tractography entities (i.e. collections of streamlines), _but which may nonetheless correspond to the same anatomical (and thus the same named) structure_ .

#### resortClassificationStruc

This function resorts an input [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) such that left and right variants of a tract are next to eachother.  Future versions may attempt to sort singleton tracts to the front and subcategories (i.e. interhemispheric or ufibers) in an additional fashion.

#### classificationStrucGrouping

Conceptually similar to [resortClassificationStruc](resortclassificationstruc), but merges them into a non-hemisphere specific category.  Uses a more birtle method than [resortClassificationStruc](resortclassificationstruc), [likely needs update](https://github.com/DanNBullock/wma_tools/issues/22)

#### removeOutliersClassification

A [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b)-specific implementation of [mbaComputeFibersOutliers](https://github.com/francopestilli/mba/blob/master/compute/mbaComputeFibersOutliers.m), a cleaning method for tracts which removes streamlines that are some number of deviations to long or too far away from the tract centroid.  Iterates across the tract identies contained within the [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) and performs the cleaning step with the same parameter settings for all.

#### clearNonvalidClassifications

Using a pairing of a [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) and a [LiFE structure](https://brainlife.io/datatype/58d15eaee13a50849b258844) this function removes those streamlines from the classification labeling schema such that only those streamlines with evidence maintain their identity.

### Debug_Tools

The functions contained within the **Debug_Tools** directory are useful for interactive segmentation and segmentation diagnostics

#### plotEndpointsOnFG

A plotting function that plots a collection of streamlines (entered as an FG) along with their endpoints and midpoints (in separate colors)

#### plotROIandFG

This function plots both a fiber group and an ROI in the same reference space.  Helpful for identifying proximity of an ROI to streamlines.

#### quickPlot

Lightweight streamline plotting function that undelies the streamline debugging visualizations.

#### quickPlotClassByName

iterates across the entries in a [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) and quickly plots the associated structures.  Good for doing QC on a segmentation.  Be wary of using this on a classification structure with a large number of categories/names.

### ROI_Tools

This directory contains the "meat" of the wma_tools library.  Here, there are the majority of the ROI extraction and modification functions that facilitate wma_tools anatomically guided segmentations.

#### amalgumROIsFromInput

This function extracts either amalgum (in that there are multiple requested labels), or singular (in that there is a single label) ROIs fom an input atlas/parcellation NIfTI; **OR** if a input directory is passed, this function will produce either singular (in that a single ROI file is requested; i.e. simple passthrough) or amalgum (in that multiple ROI files are merged) output ROI structures.  Makes no assumptions about ROI file name conventions (requires **exact** name input as string, other than ".nii.gz", which the function appends).  Features an input parameter for inflation.  Useful for creating pipelines that are interoperable with methods that have used  ADFGAG

#### MultiIntersectROIs

This function finds the overlap of multiple ROIs.  If an atlas is passed in, then the **varargin** can be used to specify which label overlaps you seek (presumably after non-trivial inflation, which is also a parameter input).  Alternatively, if **varargin** is a Series of mask ROIS, simply outputs the intesection  of those ROIs.

#### amalgamateROISacrossDirectories

This function iterates across input ROI directories and creates an amalgum across all ROIs of the same name across input directories.  IMPORTANT NOTE: given that this is a simple amalgamation of ROIs, it is presumed that all ROIs of the same name are in the same reference space and are the same dimensions.  This would be used, for example, once a multiple subjects' atlases or ROI directories had been warped to some template space, and then some amalgamtion needed to be done.  This was [implmented here](https://github.com/DanNBullock/EcogAnalysisCode) for [Babo-Rebelo et al. 2020](https://www.biorxiv.org/content/10.1101/2020.05.22.102046v1)

#### atlasROINumsFromCoords

Determines the atlas label associated with the input coordinate.  Input coordinate can be either in image space or in acpc (e.g. affined) space.

#### endpointAtlasCriteria

This function determines which streamlines end in the spcecified labels of the input atlas.  This function is not recomended for extremely large groups of streamlines.

#### intersectROIs

Finds the simple intersection of two point cloud ROIs (i.e. [vistasoft](https://github.com/vistalab/vistasoft) format).  No distance compuation performed, thus point matches must be exact.

#### loadAndParseROI

Takes an input string or object that corresponds to either a vistasoft ROI or a NIfTI ROI and outputs a standard, [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI object

#### makePlanarROI

Creates a planar ROI, in [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI format, at the specified, coordinate space (i.e. **not** image space) location PROVIDED AS A SINGLE INTEGER (i.e. -20 mm).  Input parameter "dimension" indicates the "stable" dimension of the plane, i.e. the coordinate that will be shared amongst all points associated with the planar ROI.  ACPC orentation recommended in order to ensure that planes are orthogonal to standardly oriented anatomy.

#### mergeROIs

Merges the coordinates of two [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI objects into a single ROI output.

#### midpointROISegment

This function determines whether the midpoint of a set of streamlines intersects with an ROI.  Useful for performing a midpoint specific segmentation step.  Midpoints can be input directly, or can be computed from input tractogram.  In either case output is a boolean vector indicating streamline-specific intersection with the (format agnostic, i.e. NIfTI or [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI objects) input ROI.

#### modifyROI

This function modifies an input ROI (either as direct input of a [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI objects, or as specified by an atlas-label(sequence) pairing) by removing all ROI coordinates NOT specified by the input.  For example, relative to a particular 3-D input coordinate (i.e. x,y,z), an input instruction (in the "location" variable) of "superior" will preserve all ROI cordinates superior of the specified Z value of the input coordinate, thereby removing all coordinates inferior to the specified refCoord.  Notably the input **refCoord** can actually be an planar [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI object.  The desired coordinate will be inferred from the singleton (i.e. "stable") value of the ROI.

#### planeFromROI

This function creates a creates a planar ROI, in [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI format, at the specified boundary of the input ROI.  Input ROI format can be either a direct input of a [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI object, or as specified by an atlas-label(sequence) pairing.  Desired boundary of ROI is specfied by string input from the following options ("top", "superior", "bottom", "inferior", "anterior", "posterior", "medial")

#### roiFromAtlasNums

This function extracts an ROI, in [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI format, from the input atlas.  Desired atlas labels are specified using the **ROInums** input parameter, and can be a vector of multiple integer values (which will result in a single, merged ROI).  **smoothKernel** parameter input permits inflation of extracted ROI.

#### subtractROIs

This function subtracts the first input ROI, in [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI format, from the second input ROI, in [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI format.

#### tractByEndpointROIs

This function segments streamlines from an input fg/tractogram such that, for all surviving stramlines one streamline endpoint is in each of the input rois.  Prevents odd situations that were arising with other segmentation methods (i.e. inclusion of within-roi u fibers).  Input ROIS must be [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI format and formatted thusly:  [{ROI1} {ROI2}]

#### dtiRoiFromNiftiObject
A copy of the [vistasoft](https://github.com/vistalab/vistasoft) by the same name, however this function doesn't presume that the NIfTI being passed is a string to a saved file.  This permits rapid, iterated use without massive load times.  Overall, this function converts a standard NIfTI formatted mask to a [vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI.  Used to obviate [roiFromAtlasNums](#roifromatlasnums) with [dtiRoiFromNiftiObjectSmoothWrapper](dtiroifromniftiobjectsmoothwrapper)

#### dtiRoiFromNiftiObjectSmoothWrapper

Essentially, the same function as [dtiRoiFromNiftiObject](dtiroifromniftiobject), except that it removes unneeded input and incldues a parameter for inflation in order to replicate functionality of [roiFromAtlasNums](#roifromatlasnums).

#### multiROIpairSeg

This function segments tracts from an input wbfg by teratively applying [tractByEndpointROIs](#tractbyendpointrois), which segments tracts by their endpoints.  In theory, input ROIs can be either .mat ([vistasoft](https://github.com/vistalab/vistasoft) point cloud ROI) format or NIfTI format, either as strings or as
 objects.  This function is essentially the base function for [SegROIfromPairStringList](#segroifrompairstringlist).
 
#### SegmentFascicleFromConnectome

An adapted version of [feSegmentFascicleFromConnectome](https://github.com/francopestilli/life/blob/9a169efb9219137c8b0b43fbf4261922d287fe36/utility/feSegmentFascicleFromConnectome.m) with various modifications.  Takes an input fiber group (i.e. tract, tractome, whole brain fiber group) and iteratively applies ROI based criteria (i.e. 'NOT' or 'AND') or more specific options (i.e. 'endpoints', which can have its output negated after running) to output both a new tract and a boolean vector reflecting the survivors.

### Segmentations

This directory contains various instantiations of segmentations implemented with wma_tools.  Given the iterative and archival nature of the contents of this directory, only certian notable examples will be discussed.

#### genericSegmentFunction

This function serves as the template for the new, exhaustively documented format for wma_tools based segmentations.  It divides up any single tract segmentation into 4 primary sections, some of which have their own subsections 

- 1. the establishment of a category criterion for this tract (using [streamlineCategoryPriors](#streamlinecategorypriors)), which indicates the general brain lobes involved (e.g. "posterior" , "parietal, etc.)
- 2. the establishment of additional endpoint criteria 
    - 1.  Extract the relevant rois from atlases into storage objects, and modify them if necessary
    - 2.  Segment streamline candidates using the refined endpoint ROI criteria
    - 3.  Apply additional endpoint, related criteria if necessary.
- 3. application of generic, anatomically informed criteria
    - 1.  Generation and application of anatomically informed planes
    - 2.  Generation and application of anatomically informed volumetric ROIs
    - 3.  Generation and application of other anatomically informed criteria (i.e. midpoint criteria)
- 4.  The compilation of the various criteria in to a single omnibus criteria which is used to generate this tract's [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) entry.

##### segmentAntPostTracts

The latest version of this function (v4) is a work-in-progress look at the application of the [genericSegmentFunction](#genericsegmentfunction) template.

####  segmentSuperficialFibers

This function contitutes an attempt at segmenting streamlines which remain close to the grey matter-white matter border.  These are typically associated with U-fibers.  Previous research has indicated that tracts with longer traversals are deeper in the white matter

#### streamlineCategoryPriors

This segmentation is a workhorse function which exhaustively assigns streamlines identies in accordance with their terminations.  The possible categories are 'subcorticalROIS','spineROIS','cerebellumROIS','ventricleROIS','wmROIS','ccROIS','unknownROIS','OpticCROI','FrontalROIs','TemporalROIs','OccipitalROI','ParietalROI','pericROI','insulaROI','thalamicROI','caudateNAcROI','lenticularNROI','hippAmig'.  Because all streamlines have two terminations, all streamlines are assigned a pairing of categories.  This function is used to prefilter streamlines for most wma_tools segmentations.  Additionally this function can provide a decent quality assurance check for a tractome output, in order to assess general tractome features [as described here](https://github.com/brainlife/app-streamlineCategorySegmentation/tree/2.0)

This segmentation can be run as a standalone function on brainlife.io:  
[![Run on Brainlife.io](https://img.shields.io/badge/Brainlife-bl.app.249-blue.svg)](https://doi.org/10.25663/brainlife.app.249)

### Stream_Tools

This directory contains functions used to modify, interact, or quantatively assess streamlines.  Primarily uses [vistasoft](https://github.com/vistalab/vistasoft) fg format for tractography

#### applyEndpointCriteria

This function evaluates streamlines from the input FG against a number of criteria.  These criteria are specified in a specific fashion using the **varargin** variable input.  Allows for some relatively complex and advanced requests.  Specifically:

> the varargin is interpreted in **sequential triplets** of inputs such that
-    Input (1) is a coordinate (or a planar ROI)
-    Input (2) is a relative criteria expressed as a string (i.e. "top", "superior", "bottom", "inferior", "anterior", "posterior", or "medial")
-    Input (3) is a string input of 'both', 'one', or 'neither' indicating which endpoints this criteria should be true for.

Outputs a boolean vector indicating satisfaction of all criteria.

#### applyMidpointCriteria

Similar to [applyEndpointCriteria](#applyendpointcriteria), except that **varargin** variable input is interpreted in **pairs**, as there is no need to specify which "endpoints" the criteria should be true of (there is only one midpoint).  As was the case for [midpointROISegment](#midpointroisegment) streamline-related inputs can either be in the form of pre-computed midpoints (for computational speedup) or an input FG structure.

#### classifiedStreamEndpointCortex

This is a wrapper around [endpointMapsDecay](#endpointmapsdecay) which iteratively generates (and saves) NIfTI-formatted termination maps for each structre identified in an input [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  subSelect input parameter permits sub-selection (by index) of those structure names which the user wishes to extract.  decayFunc and decayRadiusThresh parameter inputs control application of smoothing algorithms.

#### investigateTract

This function applies [atlasROINumsFromCoords](#atlasroinumsfromcoords) to both endpoints of each streamline from the input FG, and thereby determines which two labels from the input atlas are assocated with each streamline.

#### orientFgRAS

This function determines the "primary dimension of traversal" (i.e. the dimension the streamline travels the most distance in) for each streamline and reorients the streamlines in accordance with RAS-LPI conventions such that the first endpoint is the RAS endpoint and the second endpoint is the LPI endpoint.  Dimension of assesment is determined by "primary dimension of traversal".  Functionality underwritten by [endpointClusterProto](#endpointclusterproto)

#### reorientStreamline

Single-streamline based implementation of [orientFgRAS](#orientfgras)

#### reorientFiber

Predecessor function to [orientFgRAS](#orientfgras), before the implementation of [endpointClusterProto](#endpointclusterproto)

#### splitTractAtPoint

This function splits an input FG-formatted tract at a specified coordinate, based on the relative location of the closest node of a given streamline to a given
plane.  As an illustration, this code is the generalized version of previous code which was used to separate the pArc and TPC from the same body at the bottom of the IPS for versions of the segmentation associated with [Bullock et al. 2019](https://doi.org/10.1007/s00429-019-01907-8).  **coordinate** parameter input is an X Y Z coordinate, while **location** parameter input is a string ('x', 'y', or 'z') indicating the dimension for assesment.

#### endpointClusterProto

An over-engineered function which uses something akin to the [MinimumAverageDirectFlipMetric](https://dipy.org/documentation/0.16.0./reference/dipy.segment/#minimumaveragedirectflipmetric) as implemented by [Garyfallidis et al. 2012](https://dx.doi.org/10.3389%2Ffnins.2012.00175) to classify streamline endpoints as being associated with either the RAS or LPI end of the tract.  Outputs both the endpoint identities and the coordinates.

### Utils

The **Utils** directory contains utility functions that are not specific to any of the above categories.  They will not be covered invidually in this readme.

#### External

Contains functions obtained from outside sources, all rights are reserved by their respective authors

#### GitUtils

This directory contains a number of functions that were used to manage the development of the wma_tools library, specifically as it relates to github.

#### Loaders

Loaders features functions that were used to load datatypes which could be stored in myriad formats.

### Visualization

This directory contains code related to visualization

#### plotClassifiedStreamsAdaptive

This function is used to generate publication quality streamline anatomy figures from an input classification structure.  The funtion plots the anatomy against a t1 cross section (of the specified view) using tools from the [Matlab Brain Anatomy](https://github.com/francopestilli/mba) toolkit.  The **subSelect** input parameter allows users to specify the indexes of those classification structure entries that they wish to plot, while the input **colors** alllows users to specify the colors for those structures.

#### MBAplotTracts_BL

This function is a brainlife.io wrapper for [plotClassifiedStreamsAdaptive](plotclassifiedstreamsadaptive), and was previously associated with a standalone app.
