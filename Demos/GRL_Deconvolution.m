%% Perform GRL spherical deconvolution (https://arxiv.org/abs/1910.05372)
% MyMat is a .MAT file in ExploreDTI-like format and containing a multi-shell acquisition.
% A guide on how to reach the .MAT is available in CSD_Tractography_Script.m

% Converts a .MAT to the MRIToolkit format
mrt_data = EDTI.EDTI_Data_2_MRIToolkit('mat_file','mymat.mat');
% Prepare the spherical deconvolution class on this dataset
SD = SphericalDeconvolution('data',mrt_data);
% Add an anisotropic response function using the DKI model
SD.AddAnisotropicRF_DKI([1.7e-3 0.2e-3 0.2e-3],0); % WM
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
SphericalDeconvolution.SaveOutputToNii(SD,GRL_Results,'GRL_Test');

%% Perform mFOD spherical deconvolution (https://www.biorxiv.org/content/10.1101/739136v1)

mrt_data = EDTI.EDTI_Data_2_MRIToolkit('mat_file','mymat.mat');
SD = SphericalDeconvolution('data',mrt_data);
% Add an anisotropic response function using the DKI model
SD.AddAnisotropicRF_DKI([1.7e-3 0.2e-3 0.2e-3],0); % WM
% Add an anisotropic response function using the NODDI model
SD.AddAnisotropicRF_NODDI({[0.4 1.7E-9 1 0.0 3E-9 1]}); % GM
SD.AddIsotropicRF(3e-3); % CSF
SD.setInnerShellWeighting(1.0); % in mFOD this parameter has not been investigated
SD.setDeconvMethod('L2LSQ');
mFOD_Results = SD.PerformDeconv();
SphericalDeconvolution.SaveOutputToNii(SD,mFOD_Results,'mFOD_TEST');
