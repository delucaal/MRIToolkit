% Example of processing of dMRI data to perform constrained spherical deconvolution (CSD)
% and the subsequent fiber tractography. For this example I am going to assume HCP-like data format

% Just needed to compute the b-matrix
copyfile('bvals','data.bval');
copyfile('bvecs','data.bvec');

% If you want all subsequent output to be saved as compressed NIFTI
MRTQuant.EnforceNiiGz(true);
% b-matrix conversion
MRTQuant.b_Matrix_from_bval_bvec('data.bval');
% The ExploreDTI file convention disregards the Q/S matrix of the NIFTI. As a result, we actually need to adjust the data matrix
MRTQuant.FlipPermuteSpatialDimensions('nii_file','data.nii','flip',[0 1 0],'permute',[1 2 3]);
MRTQuant.FlipPermuteSpatialDimensions('nii_file','nodif_brain_mask.nii','flip',[0 1 0],'permute',[1 2 3]);
% Convert the data to .MAT format (ExploreDTI compatible). On other datasets, you should check whether the gradients flip/permute is correct. I refer you to the ExploreDTI manual for this.
MRTQuant.PerformDTI_DKIFit('nii_file','data_FP.nii','grad_perm',2,'grad_flip',2,...
    'fit_mode','ols','mask','nodif_brain_mask_FP.nii','txt_file','data.txt');

% Perform CSD to derive the Fiber Orientation Distribution (FOD)
MRTQuant.PerformCSD('mat_file','data_FP.mat','output','data_FP_CSD_FOD.nii','mscsd',0);
% Perform CSD-Based fiber tractography
MRTQuant.PerformFODBased_FiberTracking('mat_file','data_FP.mat','fod_file','data_FP_CSD_FOD.nii',...
    'SeedPointRes',[1.25 1.25 1.25],'StepSize',0.6,'FODThresh',0.3,'AngleThresh',45,...
    'FiberLengthRange',[30 500],'output','data_FP_Tracts.mat');
