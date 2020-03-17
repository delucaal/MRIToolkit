% 1 Download data from https://surfdrive.surf.nl/files/index.php/s/kAYfZ0xugHKOC4t

% Set SampleData/RawNii as working directory

% If you like .nii.gz, excute the following. 
% Warning: Always specify the input files as .nii, even when using .gz
EDTI.EnforceNiiGz(true);

% Generate a bmat .txt file from .bval - .bvec
EDTI.b_Matrix_from_bval_bvec('bval_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.bval',...
    'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt');

% Make sure the spatial orientations are consistent
% T1 first
EDTI.FlipPermuteSpatialDimensions('nii_file','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain.nii',...
    'output','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP.nii');
% Downsample the T1 at 2mm resolution to save memory (not needed if your
% workstation has 16+GB ram)
EDTI.ResampleDataSpatially('nii_file','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP.nii'...
    ,'output','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP_ds.nii',...
    'res',[2 2 2]);
% Diffusion after
EDTI.FlipPermuteSpatialDimensions('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.nii',...
    'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.nii','flip',[0 1 0]);

% % Optional pre-processing steps - Signal drift correction
% copyfile('sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt',...
%     'sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.txt');
% EDTI.PerformSignalDriftCorrection('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.nii',...
%     'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_sdc.nii');

% Optional pre-processing steps - denoising with MP-PCA
MRTD.MPPCADenoising('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP.nii',...
    'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP');

% % Optional pre-processing steps - Gibbs Ringing correction
% copyfile('sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt',...
%     'sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.txt');
% EDTI.PerformGibbsRingingCorrection('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.nii',...
%     'output','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_GR.nii');

% Convert the diffusion data to ExploreDTI-like .mat
EDTI.PerformDTI_DKIFit('nii_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.nii',...
    'txt_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1.txt');

% Perform MoCo-EPI correction and final DTI fit
EDTI.PerformMocoEPI('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised.mat',...
    'epi_tgt','sub-MRI_ses-1/anat/sub-MRI_ses-1_sT1W_3D_TFE_1_brain_FP_ds.nii','fit_mode','wls');

% Export DTI metrics to .nii
EDTI.MatMetrics2Nii('sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat');

% Convert EDTI-like data to the MRIToolkit standard structure
mrt_data = EDTI.EDTI_Data_2_MRIToolkit('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat');

% Use the Spherical deconvolution class to perform the Generalized Richardson
% Lucy deconvolution

% Prepare the spherical deconvolution class on this dataset
SD = SphericalDeconvolution('data',mrt_data);
% Add an anisotropic response function using the DKI model
SD.AddAnisotropicRF_DKI([2.1e-3 0.e-3 0.e-3],0); % WM
% Add an isotropic response function with the ADC model
SD.AddIsotropicRF(0.7e-3); % GM
% Add an isotropic response function with the ADC model
SD.AddIsotropicRF(3e-3); % CSF
SD.setInnerShellWeighting(0.2); % Do not loose angular resolution due to the lower shells
SD.AutomaticDRLDamping(); % see Dell'Acqua 2013
SD.setDeconvMethod('dRL'); % damped Richardson Lucy
% Perform the actual deconvolution
GRL_Results = SD.PerformDeconv();
% Export everything to NIFTI
SphericalDeconvolution.SaveOutputToNii(SD,GRL_Results,'sub-MRI_ses-1/dwi/GRL_deconv');

EDTI.PerformFODBased_FiberTracking('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'fod_file','sub-MRI_ses-1/dwi/GRL_deconv_CSD_FOD_scaled.nii',...
    'SeedPointRes',[2 2 2],'AngleThresh',30,'StepSize',1,...
    'output','sub-MRI_ses-1/dwi/GRL_deconv_Tracking.mat');

% Filter at the GM/WM interface
SphericalDeconvolution.TerminateTractsWithFraction('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'tract_file','sub-MRI_ses-1/dwi/GRL_deconv_Tracking.mat',...
    'mask_mode','wm','fraction_file','sub-MRI_ses-1/dwi/GRL_deconv_fractions.nii',...
    'out_file','sub-MRI_ses-1/dwi/GRL_deconv_Tracking_wmborder.mat');
% Filter at the GM/CSF interface
SphericalDeconvolution.TerminateTractsWithFraction('mat_file','sub-MRI_ses-1/dwi/sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'tract_file','sub-MRI_ses-1/dwi/GRL_deconv_Tracking.mat',...
    'mask_mode','gm','fraction_file','sub-MRI_ses-1/dwi/GRL_deconv_fractions.nii',...
    'out_file','sub-MRI_ses-1/dwi/GRL_deconv_Tracking_gmborder.mat');
