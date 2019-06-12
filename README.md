 # MRIToolkit
 
 ## What is it?
 For my [my research](https://www.isi.uu.nl/people/alberto-de-luca/) I developed a set of scripts and tools to process (diffusion) magnetic resonance imaging (MRI) data, that are (slowly) being integrated in MRIToolkit in a MATLAB (R) toolbox .

## Where do I find it?
[Here](https://github.com/delucaal/MRIToolkit) on Github!

## What can I do with it?
The toolbox is slowly but steadly growing!
###### Ready to use:
- [x] **'LesionEditor'**: a graphical user interface to visualize and edit segmentations of 3D MRI images, originally designed for delination of multiple-sclerosis lesions on fluid attenuated inversion recovery (FLAIR) images. Requires MATLAB R2018a or newer. **Documentation coming soon**
- [x] **'NiftiIO_basic'**: Basic Nifti input/output, originally written by [Jimmy Shen](https://nl.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
- [x] **'OptimizationMethods'**: Classes for numeric optimization methods
- [x] **'Relaxometry'**: Classes for T1/T2 quantification using inversion-recovery / multi-echo data
- [x] **'ThirdParty'**: Utilities from third parties used in other scripts. Includes EPG code from [Brian Hargreaves](http://web.stanford.edu/~bah/software/epg/)
###### Being integrated and coming soon:
- [ ] **'ImageRegistrations'**: Image registration utils based on Elastix
- [ ] **'DW_basic'**: Basic data structures / utilities for diffusion MRI data
- [ ] **'Diffusion_basic'**: Class for (basic) dMRI quantification
- [ ] **'DW_IVIMDTDK_I'**: Diffusion MRI fit utilities - IVIM, DT, DKI 
- [ ] **'Dicom_utils'**: Tools for handling unconventional or buggy DICOMs
- [ ] **'DW_Laplacian_NNLS'**: Tools for spectral multi-compartment fit (NNLS/L2NNLS/PL2NNLS)
- [ ] **'MixedCodeUtils'**: 'Useful general purpose functions
- [ ] **'MRIfoundation'**: Classes for MRI sequences abstraction
- [ ] **'EPG_simulator'**: Classes for EPG simulations

###### Not yet planned for release:
- **'dfMRI'**,'Diffusion fMRI utilities

###### Working examples:
- **'EPGFitT2Muscle.m'**: in Relaxometry, this function fits bi-exponential T2 (water/fat) with EPG
- **'DESPOT12Fit.m'**: an example script to fit DESPOT1 and DESPOT2 as in [Deoni et al. 2005](https://www.ncbi.nlm.nih.gov/pubmed/15690526)
- **'MultiEchoTFit.m'**: an example script to fit mono-exponential T2 with multi-echo data. More advanced fit approaches (bi-exponential, EPG) will come soon
- **'InversionRecoveryT1Fit.m'**: an example script to fit T1 with inversion recovery data.

###### Notes:
- This code is a work in progress. It will be updated without notice to ensure bug-fixes and the inclusion of best available methods
- Most code is poorly commented and not general, but will be improved over releases. 

