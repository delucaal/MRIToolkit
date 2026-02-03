New in this version [1.7]:
- Initial support for calculating structural connectivity (based on tools from EDTI)

Update 11-12-2024:
- Bug fixes 
- Updated compatibility to the latest version of WhiteMatterAnalysis
- Initial support for tractography with Scilpy tools

Update 17-02-2023:
- This is the last stable release before a major change in structure. From the next release, all main functionalities will be come available as standalone commands and be designed to fully rely on Niftis only.
- Fixes to the handling of Nifti headers
- Automated brain extraction using image registration
- added new class MRT_Library that will collect all own methods currently in MRTTrack and MRTQuant (transition in progress). This will clean MRTTrack and MRTQuant which are meant as documented interfaces.
- Dramatic speed up of GRL and mFOD (major revision of MRTTrack.PerformDeconv)
- Changes to MRTQuant.LoadNifti and SaveNifti, which now allows to preserve the original header of NIFTIs (if required). In future releases, no modifications of the headers will be required.
- Update of Neuro.m to support the latest version of Elastix
- Compressed niftis are now (un)compressed in a temporary folder
- Added functions to track multiple FODs simultaneously
- Changes to the automatic handling of Q-S form of NIFTIs
- New mFOD fiber tracking
- New TerminateTractsWithFraction approach with 1 voxel tolerance
- Moved Trackers to aanother location
- Fixes to SHPrecomp
- Fixes to MRIToolkitInit
- Tractography can now use parallel computing
- Handling of NIFTI scaling factors
- Failsafe mode for registration
- Support for FSL Eddy
- Support for FOD normalization
- One call command to run GRL
- Fix for MSCSD (probably still needs amendment on the detection of data shells)
- Fix fot CAT12 processing
- command line fixes
- standardized preproc includes conformation
- fixes to conformation
- support for offset in DKI for GRL
- fixes to DKI export


Update 09-2021:
- [NEW] Checkout the [tutorial](https://github.com/delucaal/MRIToolkit/tree/master/GettingStarted/Diffusion/7_WMA) on how to install/use automated tractography clustering (WMA) algorithm with CSD/GRL/mFOD
- Fixes to MSCSD-related functions
- Fixes to WM automated clustering
- Support for Q-form in MRTTrack.ConformSpatialDimensions() 
- New FOD scaling mode for GRL/mFOD
- Renaming of some command line functions
- Added an option to customise the number of dRL iterations
- New functions to estimate SNR from b=0s/mm2 images (MRTQuant)


Update 17-02-2021:
- Run the CAT12 automatic pipeline for T1 images
- Code re-organization into two main classes: MRTQuant (preprocessing/DTI/DKI) and MRTrack (Tractography related)
- Support for automatic fiber clustering (see below for reference)
- Early support for integration with Python (needed for the point above)
- Support for VTK poly data 4.2 and 5.1 (also needed for the clustering)
- Added a robust option to GRL/mFOD deconvolution
- Initial support for storing the NIFTI Q/S form (to improve interoperability with other tools, not implemented yet)
- Integration with CAT12