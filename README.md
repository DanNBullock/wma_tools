

# White Matter Anatomy tools (wma_tools)

## Overview description

This repository contains a number of matlab code tools which can be used to perform anatomically informed segmentations and analyze the outputs.  Additionally, it includes several examples of _actually_ implemented segmentations for established white matter tracts.

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
| Atlas_tools               | Tools for working with NIfTI atlases/parcellations; both processing and analysis                                                                  |
| BL_Wrappers               | Wrappers for interacting with brainlife.io apps                                                                                                   |
| ClassificationStruc_Tools | Tools for working with [White matter classification (WMC) structures](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b)                     |
| Debug_Tools               | Tools for troubleshooting and developing segmentations; typically lightweight streamline visualizations.                                          |
| ROI_Tools                 | Tools for obtaining, modifying, and utilizing ROIs.  Both NiFTI and [vistasoft](https://github.com/vistalab/vistasoft), point-cloud format.       |
| Segmentations             | Implemented segmentations using wma_tools methods; contains single tract and multi-tract implementations; contains archived segmentation versions |
| Stream_Tools              | Tools for assessing streamline characteristics, superficial modification, and criteria application                                                |
| Utils                     | General use utility functions, externally sourced functions                                                                                       |
| Visualization             | Functions for generating visualizations (typically publication quality)                                                                           |

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
| [efficiencyRat](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L49)             | Efficiency ratio: the ratio between the displacement and the actual traversal length of each streamline.  1 = maximal efficiency, 0 = minimal efficiency (e.g. ending up exactly where it started)                                                                                                 |
| [AsymRat](https://github.com/DanNBullock/wma_tools/blob/f989c66dbfcc0bb943e9f1ed0954a6d3b52ab32a/Analysis/ConnectomeTestQ_v2.m#L52)                   | The square of the difference between the the **efficiencyRat** for the two halves of the streamline.  1 = maximal asymmetry, 0 = minimal asymmetry. \n  _Example 1_ : if the first half was completely circuitous and the second half was a straight line;  (0-1)^2 =1 \n _Example 2_ : both halves are a straight line; (1-1)^2 =0 ;                  |
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

Generates [smoothed/inflated](https://www.mathworks.com/help/matlab/ref/smooth3.html) NIfTI **density** (i.e. count) masks **for the streamline endpoints** of tracts identified in the input [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  Requires a reference NIfTI (to establish mask NIfTI dimensions) and a refrence tract/tractome (to extract identified streamlines from).  Relies on **endpointClusterProto** in order to ensure streamline endpoints are associated with appropriate endpoint group/NIfTI output.

#### quantAllWMNorm
Relevant data domain(s):  tractography

Quantifies more common quantative traits associated with white matter.  Underwrites [feAndAFQqualityCheck](#feandafqqualitycheck).  [quantAllWMNorm](#quantallwmnorm) is predicated upon the assumption of an associated [white matter classification (WMC) structure](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b) and computes a number of stats relative to the (presumed) source whole-brain tractome.  This is as opposed to [quantWBFG](#quantwbfg) which makes no presumptions regarding associated [white matter classification (WMC) structures](https://brainlife.io/datatype/5cc1d64c44947d8aea6b2d8b).  
See [this readme](https://github.com/brainlife/app-tractographyQualityCheck/tree/1.2) for an explicit listing, description, and code linkage of each quantative feature computed.

#### quantWBFG
Relevant data domain(s):  tractography

The standalone quantification function for a singleton (presumed whole-brain tractome) tract structure.  Quantifies more common quantative traits associated with white matter.  Associated with [quantAllWMNorm](#quantallwmnorm).  
See [this readme](https://github.com/brainlife/app-tractographyQualityCheck/tree/1.2) for an explicit listing, description, and code linkage of each quantative feature computed.


