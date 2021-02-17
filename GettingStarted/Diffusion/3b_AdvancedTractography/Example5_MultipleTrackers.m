% This example shows how to perform fiber tractography with the different
% methods implemented in MRIToolkit. The following lines apply to the
% example provided data after spherical deconvolution reconstruction. 

% The file SEED_CT.nii refers to a seeding mask. Either annotate your
% favorite region of interest for seeding, or remove the parameter 'SeedMask' to perform whole brain
% tractography.

global MRIToolkit

% If the following is set to empty, the standard deterministic tractography
% implemented in ExploreDTI will be used.
MRIToolkit.fibertracker.type = '';
MRTQuant.PerformFODBased_FiberTracking('mat_file','sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'fod_file','GRL_deconv_CSD_FOD_scaled.nii','StepSize',1,'FODThresh',0.01,...
    'output','GRL_DetTrackingTest_CT.mat','FiberLengthRange',[30 500],'AngleThresh',50,...
    'SeedPointRes',[2 2 2],'SeedMask','SEED_CT.nii')

% The following tracker attempts to perform multi-level fiber tractography
% as described in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4038925/ and
% in Zhylka et al. (submitted)
MRIToolkit.fibertracker.type = 'MultiPeakTracker';
% The number of levels indicates the number of times in which deterministic
% tractography is performed from the points of an existing streamline that
% contain multiple fiber orientations. Every new level extends the previous
% streamlines with new reconstructions, but the proble explodes
% computationally for nlevels > 2.
MRIToolkit.fibertracker.parameters.nlevels = 1;
MRTQuant.PerformFODBased_FiberTracking('mat_file','sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'fod_file','GRL_deconv_CSD_FOD_scaled.nii','StepSize',1,'FODThresh',0.01,...
    'output','GRL_MultiPeakTrackingTest_CT.mat','FiberLengthRange',[30 500],'AngleThresh',50,...
    'SeedPointRes',[2 2 2],'SeedMask','SEED_CT.nii')

% The following tracker is semi-deterministic. At each iteration it can
% either follow the main orientation (classic deterministic) or sample a
% direction among the main peaks (semi-deterministic). Additionally, the step
% size and angle threshold can become stochastic by allowing for a standard
% deviation of each parameter. 
MRIToolkit.fibertracker.type = 'DistProbSHTracker';
% The number of iterations in which the tracking is repreated from the provided seed
% points
MRIToolkit.fibertracker.parameters.number_of_iterations = 100;
% The standard deviation of the step size. If 0, a deterministic step size
% is used (default)
MRIToolkit.fibertracker.parameters.stepSize_sd = 0.3;
% The standard deviation of the angle threshold. If 0, a deterministic
% angle threshold is used (default)
MRIToolkit.fibertracker.parameters.maxAngle_sd = 0.3;
% The weight mode decides whether the sampling of the peaks is performed
% uniformly (each peak has equal probability) or by weighting their
% ampltiude.
MRIToolkit.fibertracker.parameters.weight_mode = 0;

MRTQuant.PerformFODBased_FiberTracking('mat_file','sub-MRI_ses-1_dMRI_B2500_S15_MB2_v2_1_FP_denoised_MD_C_trafo.mat',...
    'fod_file','GRL_deconv_CSD_FOD_scaled.nii','StepSize',1,'FODThresh',0.01,...
    'output','GRL_DistProbTrackingTest_CT.mat','FiberLengthRange',[30 500],'AngleThresh',50,...
    'SeedPointRes',[4 4 4],'SeedMask','SEED_CT.nii')


