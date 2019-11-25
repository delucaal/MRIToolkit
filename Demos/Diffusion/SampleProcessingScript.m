% A sample dataset to try this script is available here: https://surfdrive.surf.nl/files/index.php/s/kAYfZ0xugHKOC4t
%%  DICOM to NIFTI conversion
% NiftiIO_basic.OrganizeDataBIDSLike('DICOM',pwd);

%%  Study specific definitions
root_folder = fullfile(pwd,'RawNii');
subj_name = 'sub-MRI_ses-1';
subj_source_folder = fullfile(root_folder,subj_name);
subj_dest_folder = fullfile(pwd,'derivatives',subj_name);
dmri_acquisition = {'dwi','dMRI_B2500_S15_MB2_v2'};
t1_acquisition = {'anat','sT1W_3D_TFE_1_brain'};

%% dMRI spatial flip/permute and basic pre-processing

dmri_file_list = dir(fullfile(subj_source_folder,dmri_acquisition{1},[subj_name '_' dmri_acquisition{2} '*.nii']));
if(length(dmri_file_list) ~= 1)
    error(['Ambiguity in the source folder - ' subj_source_folder ' -, for data ' dmri_acquisition{2} ', bailing out']);
end

output_basename = fullfile(subj_dest_folder,dmri_acquisition{1},dmri_file_list.name(1:end-4));
if(exist(fullfile(subj_dest_folder,dmri_acquisition{1}),'dir') < 1)
    mkdir(fullfile(subj_dest_folder,dmri_acquisition{1}));
end
EDTI.b_Matrix_from_bval_bvec(fullfile(dmri_file_list.folder,[dmri_file_list.name(1:end-4) '.bval']),[output_basename '.txt']);
copyfile(fullfile(dmri_file_list.folder,dmri_file_list.name),[output_basename '.nii']);

EDTI.FlipPermuteSpatialDimensions('nii_file',[output_basename '.nii'],'flip',[0 1 0],'output',[output_basename '_FP.nii']);
output_basename = [output_basename '_FP'];

% Optional steps - here not performed
% EDTI.PerformSignalDriftCorrection('nii_file',[output_basename '.nii'],'target_bval',1200,'target_bval_tol',100,'output',[output_basename '_sdc.nii']);
% output_basename = [output_basename '_sdc'];
% 
% EDTI.PerformGibbsRingingCorrection('nii_file',[output_basename '.nii'],'output',[output_basename '_GR.nii']);
% output_basename = [output_basename '_GR'];

EDTI.PerformDTI_DKIFit('nii_file',[output_basename '.nii'],'grad_perm',2,'grad_flip',2);
EDTI.AverageNormalizedResiduals('mat_file',[output_basename '.mat'],'normalize',1,'output',[output_basename '_residuals.nii']);

%% T1 - FLAIR spatial flip/permute

t1_file_list = dir(fullfile(subj_source_folder,t1_acquisition{1},[subj_name '_' t1_acquisition{2} '*.nii']));
if(length(t1_file_list) ~= 1)
    error(['Ambiguity in the source folder - ' subj_source_folder ' -, for data ' t1_acquisition{2} ', bailing out']);
end

if(exist(fullfile(subj_dest_folder,t1_acquisition{1}),'dir') < 1)
    mkdir(fullfile(subj_dest_folder,t1_acquisition{1}));
end

t1_basename = fullfile(subj_dest_folder,t1_acquisition{1},'T1');
EDTI.FlipPermuteSpatialDimensions('nii_file',fullfile(t1_file_list.folder,t1_file_list.name),'permute',[1 2 3],'flip',[0 0 0],'output',[t1_basename '_FP.nii']);
t1_basename = [t1_basename '_FP'];

%% dMRI MoCo and  EPI correction

EDTI.PerformMocoEPI('mat_file',[output_basename '.mat'],'epi_tgt',[t1_basename '.nii'],'constraint_epi',[1 1 1],'epi_reg_mode','fa','use_normcorr',1);
output_basename = [output_basename '_MD_C_trafo'];
% EDTI.AverageNormalizedResiduals('mat_file',[output_basename '.mat'],'normalize',1,'output',[output_basename '_residuals.nii']);
EDTI.MatMetrics2Nii([output_basename '.mat']);
%% Perform DTI-based fiber tractography

EDTI.PerformDTIBased_FiberTracking('mat_file',[output_basename '.mat'],'output',[output_basename '_Tracts_DTI.mat']);

%% Perform CSD-based fiber tractography

EDTI.PerformCSD('mat_file',[output_basename '.mat'],'output',[output_basename '_CSD_FOD.nii'],'Lmax',6);
EDTI.PerformFODBased_FiberTracking('mat_file',[output_basename '.mat'],'fod_file',[output_basename '_CSD_FOD.nii'],'output',[output_basename '_Tracts_CSD.mat']);

%% Perform GRL fiber tractography

data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',[output_basename '.mat']);
SD = SphericalDeconvolution('data',data);

% Add an anisotropic response function using the DKI model
SD.AddAnisotropicRF_DKI([1.7e-3 0.2e-3 0.2e-3],0); % WM
% Add an isotropic response function with the ADC model
SD.AddIsotropicRF(0.7e-3); % GM
% Add an isotropic response function with the ADC model
SD.AddIsotropicRF(3e-3); % CSF
SD.AddIsotropicRF(20e-3); % IVIM
SD.setInnerShellWeighting(0.2); % Do not loose angular resolution due to the lower shells
SD.AutomaticDRLDamping(); % see Dell'Acqua 2013
SD.setDeconvMethod('dRL'); % damped Richardson Lucy
% Perform the actual deconvolution
GRL_Results = SD.PerformDeconv();
% Export everything to NIFTI
SphericalDeconvolution.SaveOutputToNii(SD,GRL_Results,'GRL_Test');
