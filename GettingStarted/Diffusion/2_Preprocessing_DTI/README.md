<p align="center">
<a href="https://github.com/delucaal/MRIToolkit"> 
<img src="../../../img/MRIToolkitLogo.png" height="150"/> 
 </a> 
 </p>

# MRIToolkit - Preprocessing and DTI / DKI fit [update 14-03-2020] 
The script [Example1_Preprocessing.m](Example1_Preprocessing.m) shows how to apply some common pre-processing steps to diffusion MRI data in the .nii/.bval/.bvec format. 

- 1 To start, download [this dataset](https://surfdrive.surf.nl/files/index.php/s/kAYfZ0xugHKOC4t)
- 2 Set the working directory to the RawNii folder
- 3 To save space, it might be convenient to work with nii.gz in place of .nii. Set this preference at the beginning of your script:
```matlab
EDTI.EnforceNiiGz(true);
```
- 4 Generate a .txt b-matrix file. This is required by the ExploreDTI engine (which carries out many of these calls):
```matlab
EDTI.b_Matrix_from_bval_bvec('bval_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.bval',...
    'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt');
```
- 5 We will use the T1-weighted image to perform the EPI correction using non-linear registration. To this end, all data must be in our coordinate system. Ensure it as follows:
```matlab
EDTI.FlipPermuteSpatialDimensions('nii_file','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain.nii',...
    'output','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP.nii');
```
This step touches the .nii file discarding the orientation matrix and considering only the data matrix. Please inspect your data "_FP" and compare it to the original in a viewer such as MRIcron to check whether flips / permute are needed to achieve the proper orientation.
- 6 Resample the T1 to 2x2x2mm3 - to save memory and computation time in the registration / processing of the registered diffusion data
```matlab
EDTI.ResampleDataSpatially('nii_file','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP.nii'...
    ,'output','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP_ds.nii',...
    'res',[2 2 2]);
```
- 7 Repeat step 5 on the dMRI data. This dataset needs LR flip to be anatomically correct:
```matlab
EDTI.FlipPermuteSpatialDimensions('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.nii',
    'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.nii','flip',[0 1 0]);
```
- 8 optional step: Signal drift correction (not possible on this dataset)
```matlab
% copyfile('sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt',...
%     'sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.txt');
% EDTI.PerformSignalDriftCorrection('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.nii',...
%     'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_sdc.nii');
```
- 9 MP-PCA denoising (optional):
```matlab
MRTD.MPPCADenoising('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.nii',...
    'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP');
```
- 10 Gibbs Ringing correction (optional):
```matlab
% copyfile('sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt',...
%     'sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.txt');
% EDTI.PerformGibbsRingingCorrection('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.nii',...
%     'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_GR.nii');
```
- 11 Convert the data to the ExploreDTI-like MAT format. This also performs the DTI/DKI fit 
```matlab
EDTI.PerformDTI_DKIFit('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.nii',...
    'txt_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt');
```
- 12 Perform motion / Eddy currents correction + EPI correction by targeting the down-sampled T1:
```matlab
EDTI.PerformMocoEPI('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.mat',...
    'epi_tgt','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP_ds.nii','fit_mode','wls');
```
- 13 Save the DTI metrics to .nii:
```matlab
EDTI.MatMetrics2Nii('sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat');
```
