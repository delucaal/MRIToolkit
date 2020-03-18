<p align="center">
<a href="https://github.com/delucaal/MRIToolkit"> 
<img src="../../../img/MRIToolkitLogo.png" height="150"/> 
 </a> 
 </p>

# MRIToolkit - Spherical deconvolution and Fiber tractography [update 14-03-2020] 
MRIToolkit uses the ExploreDTI engine to perform deterministic fiber tractography, either using DTI or Spherical Deconvolution techniques. In the latter case, the fiber orientation distributions (FOD) can also be exported to the MRtrix3 or Dipy format to use their tractography algorithms.

The following examples assume that the data has been [Preprocessed](../2_Preprocessing_DTI).

- To perform DTI-based fiber tractography, type the following:
```matlab
EDTI.PerformDTIBased_FiberTracking('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'FAThresh',0.2,...
    'SeedPointRes',[2 2 2],'AngleThresh',30,'StepSize',1,...
    'output','sub-MRI_ses-1/dwi/DTI_Tracking.mat')
```

- To perform spherical deconvolution, the syntax is fairly similar to the above, but an FOD must be generated first. MRIToolkit implements CSD and MSCSD (from ExploreDTI), as well as the damped Richardson Lucy (dRL), Generalized Richardson Lucy (GRL) and mFOD algorithms.
    - To perform mFOD/GRL/dRL, follow [Example2_mFOD_GRL.m](Example2_mFOD_GRL.m)
    - To perform CSD, follow [Example3_CSD.m](Example3_CSD.m)

This will generate a .mat file in ExploreDTI-like format. You can **visualize the result** with [ExploreDTI](http://www.exploredti.com) itself, my own **Visor** app (soon available), or convert it to Trackvis / MRtrix formats (tck/trk).
