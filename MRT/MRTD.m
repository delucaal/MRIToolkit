% This class implements methods for the MRIToolkit diffusion pipeline
% Differently from other classes (as EDTI), the input / output of this
% class are Niftis and other exchangeble formats (in addition to
% ExploreDTI-compatible .mat files)
classdef MRTD < handle
    methods(Static)
        
        % Create a .MAT from .nii and .bval/.bvec for internal processing
        % in a temporary directory. Input arguments:
        % nii_file: the target .nii file
        % bval_file: the compation .bval fie
        % bvec_file: the companion .bvec file
        % grad_perm: optional, how to permute the diffusion gradients
        % grad_flip: optional, how to flip the sign of the gradients
        function temp_mat_location = QuickNiiBvalBvecToMat(varargin)
            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Need to specify the target file');
            end
            
            bval_file = GiveValueForName(coptions,'bval_file');
            if(isempty(bval_file))
                bval_file = [file_in(1:end-4) '.bval'];
            end
            
            bvec_file = GiveValueForName(coptions,'bvec_file');
            if(isempty(bvec_file))
                bvec_file = [file_in(1:end-4) '.bvec'];
            end
            
            option_value = GiveValueForName(coptions,'grad_perm');
            if(~isempty(option_value))
                perm = option_value;
            else
                perm = [];
            end
            option_value = GiveValueForName(coptions,'grad_flip');
            if(~isempty(option_value))
                flip = option_value;
            else
                flip = [];
            end
            
            while(true)
                dest_basename = fullfile(tempdir,['mrtd_' num2str(randi(500000))]);
                if(isempty(dir([dest_basename '*'])))
                    break
                end
            end
            
            EDTI.b_Matrix_from_bval_bvec('bval_file',bval_file,'bvec_file',bvec_file,...
                'output',[dest_basename '.txt']);
            
            if(~isempty(perm) && ~isempty(flip))
                EDTI.PerformDTI_DKIFit('nii_file',file_in,'txt_file',[dest_basename '.txt'],...
                    'output',[dest_basename '.mat'],'grad_perm',perm,'grad_flip',flip);
            else
                EDTI.PerformDTI_DKIFit('nii_file',file_in,'txt_file',[dest_basename '.txt'],...
                    'output',[dest_basename '.mat']);
            end
            
            temp_mat_location = [dest_basename '.mat'];
            
        end
                        
        function [DataRecon, noise_map, ncomponents] = MPPCADenoising(varargin)
            % Based on MRTRIX3 implementation of J. Veraart work on MP-PCA denoising
            % (PMID 27523449)
            % Expects data_name without extension (adds .nii)
            % Extended with Martijn's (Froeling) suggestion to weight multiple windows
            % Input arguments:
            % nii_file: the target .nii file
            % output: the output name prefix. Will add _denoised.nii
            
            coptions = varargin;
            data_name = GiveValueForName(coptions,'nii_file');
            if(isempty(data_name))
                error('Need to specify the target file (nii_file)');
            end
            save_prefix = GiveValueForName(coptions,'output');
            if(isempty(save_prefix))
                error('Need to specify the output prefix');
            end
            weight_pca = 0;
            
            tic
            vol = EDTI.LoadNifti(data_name);
            % vol.img = single(vol.img)*vol.hdr.dime.scl_slope;
            
            default_extent = 5; % size of the moving window. 5 is the default for the authors
            if(nargin < 3)
                weight_pca = 0;
            end
            
            ms2 = (default_extent-1)/2;
            [sx,sy,sz,m] = size(vol.img);
            n = default_extent*default_extent*default_extent;
            msk = find(vol.img(:,:,:,1));
            
            sigmas = NaN(length(msk),1);
            ncomponents = NaN(length(msk),1);
            DataRecon = zeros(sx,sy,sz,m);
            Visits = zeros(sx,sy,sz);
            Sigma = zeros(sx,sy,sz);
            Weight = zeros(sx,sy,sz);
            
            % Compute a value of sigma for each voxel
            for i=1:length(msk)
                [x,y,z] = points3dfromlind([sx sy sz],msk(i));
                if(x-ms2 < 1 || y-ms2 < 1 || z-ms2 < 1)
                    continue
                end
                if(x+ms2 > sx || y+ms2 > sy || z+ms2 > sz)
                    continue
                end
                xrange = x-ms2:x+ms2;
                yrange = y-ms2:y+ms2;
                zrange = z-ms2:z+ms2;
                
                local_data = double(reshape(vol.img(xrange,yrange,zrange,:),n,m));
                mean_data = mean(local_data(:));
                local_data = local_data - mean_data; % FIX
                
                [sigmas(i), ncomponents(i), U, S, V] = estimate_pca_demeaned_signal(local_data,m,n);
                
                reconstructed_signal = reconstruct_truncated_pca_signal(U,S,V,sigmas(i),m,n);
                if(reconstructed_signal == -1)
                    reconstructed_signal = local_data;
                end
                rebuilt_data = reshape(reconstructed_signal,default_extent,default_extent,default_extent,m);
                rebuilt_data = rebuilt_data + mean_data;
                
                for xi=1:length(xrange)
                    for yi=1:length(yrange)
                        for zi=1:length(zrange)
                            Visits(xrange(xi),yrange(yi),zrange(zi)) = Visits(xrange(xi),yrange(yi),zrange(zi)) + 1;
                            if(weight_pca == 0)
                                Weight(xrange(xi),yrange(yi),zrange(zi)) = Weight(xrange(xi),yrange(yi),zrange(zi)) + 1;
                                DataRecon(xrange(xi),yrange(yi),zrange(zi),:) = DataRecon(xrange(xi),yrange(yi),zrange(zi),:) + rebuilt_data(xi,yi,zi,:);
                                Sigma(xrange(xi),yrange(yi),zrange(zi)) = Sigma(xrange(xi),yrange(yi),zrange(zi)) + sigmas(i);
                            else
                                Weight(xrange(xi),yrange(yi),zrange(zi)) = Weight(xrange(xi),yrange(yi),zrange(zi)) + ncomponents(i)/r;
                                DataRecon(xrange(xi),yrange(yi),zrange(zi),:) = DataRecon(xrange(xi),yrange(yi),zrange(zi),:) + ncomponents(i)/r*rebuilt_data(xi,yi,zi,:);
                                Sigma(xrange(xi),yrange(yi),zrange(zi)) = Sigma(xrange(xi),yrange(yi),zrange(zi)) + ncomponents(i)/r*sigmas(i);
                            end
                        end
                    end
                end
            end
            
            for x=1:sx
                for y=1:sy
                    for z=1:sz
                        if(Visits(x,y,z) == 0)
                            continue
                        end
                        DataRecon(x,y,z,:) = DataRecon(x,y,z,:) / Weight(x,y,z);
                        Sigma(x,y,z) = Sigma(x,y,z) / Weight(x,y,z);
                    end
                end
            end
            
            % noise_map = zeros(size(vol.img(:,:,:,1)));
            noise_map = Sigma;
            
            components = zeros(size(vol.img(:,:,:,1)));
            components(msk) = ncomponents;
            
            t = toc;
            
            if(nargin > 1)
                OUT.VD = vol.VD;
                OUT.img = DataRecon;                
                EDTI.WriteNifti(OUT,[save_prefix '_denoised.nii']);
                OUT.img = noise_map;
                EDTI.WriteNifti(OUT,[save_prefix '_noisemap.nii']);
%                 E_DTI_write_nifti_file(components,VD,[save_prefix '_ncomponents.nii']);
            end
            
            disp(['Elapsed time ' num2str(t)]);
            
        end
        
        function PerformSpectralDeconvolution(varargin)
            % Perform the Laplacian - Spectral - Deconvolution method
            % published in De Luca et al. 2018 (https://www.ncbi.nlm.nih.gov/pubmed/30052293)
            % The input MUST have multiple diffusion weightings (4+). 
            % Input arguments:
            % nii_file: the target .nii file
            % bval_file: the corresponding .bval file
            % bvec_file: the corresponding .bvec file
            % mask_file: a mask file defining the volume to fit
            % min_bval: Minimum b-value to use
            % max_bval: Maximum b-value to use
            % Lambda: L2 value of the PL2NNLS method (0.05 by default)
            % UseLP: Use the "low_res" prior in PL2NNLS 0 or 1 (default 1)
            % ExtraSens: Tries to increase sensitivity (requires high SNR
            % data) 1 or 0 (default 0)
            % OutliersSD: Threshold for outlier rejection (default 3)
            % output: the output prefix name (suffixes and .nii will be
            % added automatically)
            if(isempty(varargin))
                help('MRTD.PerformSpectralDeconvolution');
                return
            end
            coptions = varargin;
            data_name = GiveValueForName(coptions,'nii_file');
            if(isempty(data_name))
                error('Need to specify the target file (nii_file)');
            end
            parameters.data_name = data_name;
            
            bvec_name = GiveValueForName(coptions,'bvec_file');
            if(isempty(bvec_name))
                error('Need to specify the bvec file (bvec_file)');
            end
            parameters.bvec_name = bvec_name;
            
            bval_name = GiveValueForName(coptions,'bval_file');
            if(isempty(bval_name))
                error('Need to specify the bval file (bval_file)');
            end
            parameters.bval_name = bval_name;
           
            mask_name = GiveValueForName(coptions,'mask_file');
            if(isempty(mask_name))
                error('Need to specify the mask file (mask_file)');
            end
            parameters.mask_name = mask_name;
            
            output_prefix = GiveValueForName(coptions,'output');
            if(isempty(output_prefix))
                error('Need to specify the output prefix (output)');
            end
            parameters.output_prefix = output_prefix;
            
            parameters.min_bval = GiveValueForName(coptions,'min_bval');
            if(isempty(parameters.min_bval))
                error('Need to specify the minimum b-value (min_bval)');
            end
            parameters.min_bval = parameters.min_bval;
            
            parameters.max_bval = GiveValueForName(coptions,'max_bval');
            if(isempty(parameters.max_bval))
                error('Need to specify the maximum b-value (max_bval)');
            end
            parameters.max_bval = parameters.max_bval;
            
            LambdaGrid = GiveValueForName(coptions,'Lambda');
            if(isempty(LambdaGrid))
                LambdaGrid = 0.05;
            end
            parameters.LambdaGrid = LambdaGrid;

            OutliersSD = GiveValueForName(coptions,'OutliersSD');
            if(isempty(OutliersSD))
                OutliersSD = 3;
            end
            parameters.outliers_sigma = OutliersSD;

            ExtraSens = GiveValueForName(coptions,'ExtraSens');
            if(isempty(ExtraSens))
                ExtraSens = 0;
            end
            parameters.ExtraSens = ExtraSens;

            UseLP = GiveValueForName(coptions,'UseLP');
            if(isempty(UseLP))
                UseLP = 1;
            end
            parameters.use_low_prior = UseLP;

            GMMBatch_NMRBiomed_fitscript_clean(parameters)

        end
        
        function PerformSegmentedIVIMFit2comp(varargin)
            % Perform an IVIM "segmented" fit with 2 compartments (tissue + blood).
            % The input MUST have multiple diffusion weightings (2+) including low b-values (b<=200s/mm2). 
            % Input arguments:
            % nii_file: the target .nii file
            % bval_file: the corresponding .bval file
            % bvec_file: the corresponding .bvec file
            % txt_file: the target .txt file (overrides
            % bval_file,bvec_file)
            % mat_file: the target .mat file (overrides nii_file,
            % bval_file, bvec_file)
            % mask_file: a mask file defining the volume to fit
            % min_bval: Minimum b-value to use
            % max_bval: Maximum b-value to use
            % b_thresh: The b-value after which IVIM is "negligible"
            % (typically b=200s/mm2)
            % output: the output prefix name (suffixes and .nii will be
            % added automatically)   
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = EDTI.LoadNifti(data_name);

                txt_name = GiveValueForName(coptions,'txt_file');

                if(isempty(txt_name))
                    bvec_name = GiveValueForName(coptions,'bvec_file');
                    if(isempty(bvec_name))
                        error('Need to specify the bvec file (bvec_file)');
                    end

                    bval_name = GiveValueForName(coptions,'bval_file');
                    if(isempty(bval_name))
                        error('Need to specify the bval file (bval_file)');
                    end
                    data.bvals = load(bval_name)';
                    data.bvecs = load(bvec_name)';
                else
                    gmat = load(txt_name);
                    [data.bvals,data.bvecs] = EDTI.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
            end
            
            mask_name = GiveValueForName(coptions,'mask_file');
            if(isempty(mask_name))
                error('Need to specify the mask file (mask_file)');
            end
            
            output_prefix = GiveValueForName(coptions,'output');
            if(isempty(output_prefix))
                error('Need to specify the output prefix (output)');
            end
            
            parameters.min_bval = GiveValueForName(coptions,'min_bval');
            if(isempty(parameters.min_bval))
                error('Need to specify the minimum b-value (min_bval)');
            end
            
            parameters.max_bval = GiveValueForName(coptions,'max_bval');
            if(isempty(parameters.max_bval))
                error('Need to specify the maximum b-value (max_bval)');
            end

            parameters.bthresh = GiveValueForName(coptions,'b_thresh');
            if(isempty(parameters.bthresh))
                disp('Using the default threshold b-value (bthresh)');
                parameters.bthresh = 200;
            end

            mask = EDTI.LoadNifti(mask_name);
            data.mask = mask.img;
            
            [Dhigh,Dlow,f] = SegmentedIVIMFit2Comp(data,parameters);

            EDTI.WriteNifti(f,[output_prefix '_fivim.nii']);
            EDTI.WriteNifti(Dlow,[output_prefix '_Dstar.nii']);
            EDTI.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
        end
        
        function PerformFW_Fit(varargin)
            % Perform a FW + ADC fit (ISO) with 2 compartments (tissue + blood).
            % This is a non-linear fit (slow)
            % The input MUST have multiple diffusion weightings (2+) including low b-values (b<=1000s/mm2). 
            % Input arguments:
            % nii_file: the target .nii file
            % bval_file: the corresponding .bval file
            % bvec_file: the corresponding .bvec file
            % txt_file: the target .txt file (overrides
            % bval_file,bvec_file)
            % mat_file: the target .mat file (overrides nii_file,
            % bval_file, bvec_file)
            % mask_file: a mask file defining the volume to fit
            % min_bval: Minimum b-value to use
            % max_bval: Maximum b-value to use
            % output: the output prefix name (suffixes and .nii will be
            % added automatically)   
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = EDTI.LoadNifti(data_name);

                txt_name = GiveValueForName(coptions,'txt_file');

                if(isempty(txt_name))
                    bvec_name = GiveValueForName(coptions,'bvec_file');
                    if(isempty(bvec_name))
                        error('Need to specify the bvec file (bvec_file)');
                    end

                    bval_name = GiveValueForName(coptions,'bval_file');
                    if(isempty(bval_name))
                        error('Need to specify the bval file (bval_file)');
                    end
                    data.bvals = load(bval_name)';
                    data.bvecs = load(bvec_name)';
                else
                    gmat = load(txt_name);
                    [data.bvals,data.bvecs] = EDTI.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
            end
            
            mask_name = GiveValueForName(coptions,'mask_file');
            if(isempty(mask_name))
                error('Need to specify the mask file (mask_file)');
            end
            
            output_prefix = GiveValueForName(coptions,'output');
            if(isempty(output_prefix))
                error('Need to specify the output prefix (output)');
            end
            
            parameters.min_bval = GiveValueForName(coptions,'min_bval');
            if(isempty(parameters.min_bval))
                error('Need to specify the minimum b-value (min_bval)');
            end
            
            parameters.max_bval = GiveValueForName(coptions,'max_bval');
            if(isempty(parameters.max_bval))
                error('Need to specify the maximum b-value (max_bval)');
            end

            mask = EDTI.LoadNifti(mask_name);
            data.mask = mask.img;
            
            [Dhigh,f] = FWFit2Comp(data,parameters);

            EDTI.WriteNifti(f,[output_prefix '_ffw.nii']);
            EDTI.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
            
        end                
 
        function PerformIVIM_FW_Fit(varargin)
            % Perform an IVIM + FW + ADC non-linear fit.
            % The input MUST have multiple diffusion weightings (2+) including low b-values (b<=200s/mm2)
            % as well as intermediate b-values (200 < b < 1000s/mm2). 
            % Input arguments:
            % nii_file: the target .nii file
            % bval_file: the corresponding .bval file
            % bvec_file: the corresponding .bvec file
            % txt_file: the target .txt file (overrides
            % bval_file,bvec_file)
            % mat_file: the target .mat file (overrides nii_file,
            % bval_file, bvec_file)
            % mask_file: a mask file defining the volume to fit
            % parameters.min_bval: Minimum b-value to use
            % parameters.max_bval: Maximum b-value to use
            % output: the output prefix name (suffixes and .nii will be
            % added automatically)   
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = EDTI.LoadNifti(data_name);

                txt_name = GiveValueForName(coptions,'txt_file');

                if(isempty(txt_name))
                    bvec_name = GiveValueForName(coptions,'bvec_file');
                    if(isempty(bvec_name))
                        error('Need to specify the bvec file (bvec_file)');
                    end

                    bval_name = GiveValueForName(coptions,'bval_file');
                    if(isempty(bval_name))
                        error('Need to specify the bval file (bval_file)');
                    end
                    data.bvals = load(bval_name)';
                    data.bvecs = load(bvec_name)';
                else
                    gmat = load(txt_name);
                    [data.bvals,data.bvecs] = EDTI.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = EDTI.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
            end
            
            mask_name = GiveValueForName(coptions,'mask_file');
            if(isempty(mask_name))
                error('Need to specify the mask file (mask_file)');
            end
            
            output_prefix = GiveValueForName(coptions,'output');
            if(isempty(output_prefix))
                error('Need to specify the output prefix (output)');
            end
            
            parameters.min_bval = GiveValueForName(coptions,'min_bval');
            if(isempty(parameters.min_bval))
                error('Need to specify the minimum b-value (min_bval)');
            end
            
            parameters.max_bval = GiveValueForName(coptions,'max_bval');
            if(isempty(parameters.max_bval))
                error('Need to specify the maximum b-value (max_bval)');
            end
            
            mask = EDTI.LoadNifti(mask_name);
            data.mask = mask.img;
            
            [Dhigh,Dlow,ffw,fivim] = IVIMFWFit3Comp(data,parameters);

            EDTI.WriteNifti(ffw,[output_prefix '_ffw.nii']);
            EDTI.WriteNifti(fivim,[output_prefix '_fivim.nii']);
            EDTI.WriteNifti(Dlow,[output_prefix '_Dstar.nii']);
            EDTI.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
        end        
        
    end
end

% Helper: finds a parameter by name when using varargin
function value = GiveValueForName(coptions,name)
value = [];
for ij=1:2:length(coptions)
    if(strcmpi(coptions{ij},name))
        value = coptions{ij+1};
        return
    end
end
end

% Helper function for MP-PCA denoising
function reconstructed_signal = reconstruct_truncated_pca_signal(U,S,V,sigma2,m,n)
eigv = flip(diag(S)).^2;
r = min(n,m);

sigsq1 = sigma2^2;

lam_r = eigv(1) / n;
clam = 0;
cutoff_p2 = 0;
for p=1:r
    lam = eigv(p) / n;
    clam = clam+lam;
    gam = (m-r+p)/n;
    sigsq2 = (lam-lam_r)/4/sqrt(gam);
    if(sigsq2 < sigsq1)
        cutoff_p2 = p+1;
    end
end
cutoff_p2 = r-cutoff_p2;
eigv = flip(eigv);
if(cutoff_p2 > 1)
    Snew = zeros(size(S));
    Sd = diag(sqrt(eigv(1:cutoff_p2)));
    Snew(1:cutoff_p2,1:cutoff_p2) = Sd;
else
    reconstructed_signal = -1;
    return
end

reconstructed_signal = U*Snew*V';
end

% Helper function for MPPCA denoising
function [sigma, ncomponents, U, S, V] = estimate_pca_demeaned_signal(local_data, m, n, force_sigma)
if(nargin < 4)
    force_sigma = 0;
end
r = min(n,m);
try
    [U,S,V] = svd(local_data,'econ');
catch
    sigma = -1;
    ncomponents = -1;
    U = 0;
    V = 0;
    S = 0;
    return
end

eigv = flip((diag(S)).^2);

lam_r = eigv(1) / n;
clam = 0;
sigma2 = NaN;

cutoff_p = 0;
for p=1:r
    lam = eigv(p) / n;
    clam = clam+lam;
    gam = (m-r+p)/n;
    if(force_sigma < 1)
        sigsq1 = clam/(p) / max(gam,1);
    else
        sigsq1 = sigma2^2;
    end
    sigsq2 = (lam-lam_r)/4/sqrt(gam);
    if(sigsq2 < sigsq1)
        sigma2 = sqrt(sigsq1);
        cutoff_p = p+1;
    end
end
cutoff_p = r-cutoff_p;

ncomponents = cutoff_p;

sigma = sigma2;

end

% Helper function
function [x,y,z] = points3dfromlind(MS,lin_ind)
% Linear Index For Matrix Of Size MS(3 dimensional) @ Position (x,y,z)
    nR = MS(1);
    nC = MS(2);
    z = fix((lin_ind-1)/((nR)*(nC)))+1;
    lin_ind = lin_ind-(z-1)*(nR)*(nC);
    y = fix((lin_ind-1)/nR)+1;
    x = lin_ind-(y-1)*(nR);
end

% This code performs the Laplacian - Spectral - Deconvolution method
% published in De Luca et al. 2018 (https://www.ncbi.nlm.nih.gov/pubmed/30052293)
function GMMBatch_NMRBiomed_fitscript_clean(parameters)
% Initial setup - logistics
% just name change of some parameters for convenience

data_name = parameters.data_name;
bvec_name = parameters.bvec_name;
bval_name = parameters.bval_name;

mask_name = parameters.mask_name;

parameters.min_bval = parameters.min_bval;
parameters.max_bval = parameters.max_bval;

% Load the data

data = DW_LoadData(data_name,bvec_name,bval_name,mask_name,0,'',1);
data_name_parts = strsplit(data_name,'.');
for k=2:length(data_name_parts)-1
   data_name_parts{1} = [data_name_parts{1} '.' data_name_parts{k}]; 
end

% Historically I had some issues with the b=0s/mm2 (the first volume) in my
% data, and I always used b=1 or b=5s/mm2 as starting point. Maybe the following can
% be safely commented.

data.img = data.img(:,:,:,2:end);
data.bvals = data.bvals(2:end);
data.bvecs = data.bvecs(2:end,:);

% Geometric averaging of the data per b-value

avg_data.hdr = data.hdr;
avg_data.untouch = 1;
avg_data.bvals = unique(round(data.bvals));
avg_data.img = zeros([size(data.img(:,:,:,1)) length(avg_data.bvals)]);
avg_data.bvecs = zeros(length(avg_data.bvals),3);
avg_data.mask = data.mask;

for ij=1:length(avg_data.bvals)
	SEL = abs(data.bvals-avg_data.bvals(ij)) < 1;
    avg_data.img(:,:,:,ij) = geomean(squeeze(data.img(:,:,:,SEL)),4);
end

data = avg_data;
clear avg_data;

% Filter according to initial requirements
GP = data.bvals >= parameters.min_bval & data.bvals <= parameters.max_bval;
data.img = data.img(:,:,:,GP);
data.bvals = data.bvals(GP);
data.bvecs = data.bvecs(GP,:);

% Basic definitions

params.bvals = data.bvals;
params.bvecs = data.bvecs;

% BUILD Signals Dictionary

MinD = 0.1e-3; % Minimum diffusion coefficient in solution space
MaxD = 1000e-3; % Maximum

NSignals = 300; % Number of possible solutions
DRange = logspace(log10(MinD),log10(MaxD),NSignals);

Dictionary = DW_BuildIsotropicDictionary(data,DRange);
[GMMDictionary,GMMCouples] = DW_BuildGMMDictionaryWithCR(DRange,[5 10 20 30]);

% Perform the two step deconvolution

LSQNONNEG_AVG_PROFILE = zeros([size(data.img(:,:,:,1)) length(DRange)]);
L2LSQNONNEG_PRIOR_AVG_PROFILE = zeros([size(data.img(:,:,:,1)) length(DRange)]);
LSQNONNEG_CL_AVG_PROFILE = zeros([size(data.img(:,:,:,1)) length(DRange)]);

OUTLIERS = ones(size(data.img));

op_lp = optimset('TolX',1e-2);
op_hp = optimset('TolX',1e-8);
if(parameters.ExtraSens == 1)
    op_lp = optimset('TolX',1e-4);
    op_hp = optimset('TolX',1e-16);    
end
op_e16 = optimset('TolX',1e-16);

% Step 1 - perform voxel-wise NNLS to learn the prior (with outlier
% rejection)

min_b = min(data.bvals);

disp('Warning - tolerances have been changed');

for x=1:size(data.img,1)
disp(['Total: ' num2str(x/size(data.img,1)*100) '%']);
for y=1:size(data.img,2)
for z=1:size(data.img,3)
    if(data.mask(x,y,z) == 0 || isnan(data.mask(x,y,z)))
        continue
    end
    NoisySimulatedSignal = double(squeeze(data.img(x,y,z,:)));
    if(sum(NoisySimulatedSignal == 0) > 1)
        continue
    end
    NoisySimulatedSignal = NoisySimulatedSignal/mean(NoisySimulatedSignal(data.bvals <= min_b));
    
    % This solution is more sensitive to small components but also to noise 
    [LSQNONNEG_CL_DECONV,goodpoints] = DW_RobustDeconvFinal(Dictionary,NoisySimulatedSignal,parameters.outliers_sigma,op_hp);
    OUTLIERS(x,y,z,:) = goodpoints;
    
    LSQNONNEG_CL_AVG_PROFILE(x,y,z,:) = LSQNONNEG_CL_DECONV;
   
    % This solution is much more robust to noise, and used to identify the
    % number of classes to sum up the fractions
    LSQNONNEG_AVG_PROFILE(x,y,z,:) = DW_RobustDeconvFinal(Dictionary,NoisySimulatedSignal,3,op_lp);
end
end
end

% This is the voxel-wise average of the high sensitivity solution, which
% becomes the prior
PRIOR = squeeze(nanmean(nanmean(nanmean(LSQNONNEG_CL_AVG_PROFILE(:,:,:,:))))); 
PRIOR(DRange > 500e-3) = 0;
PRIOR = PRIOR/sum(PRIOR);

% This is the voxel-wise average of the low sensitivity solution. This
% might also represent an excellent prior if you are interested only in
% signal fractions. If you would like to investigate the diffusion
% coefficient of the components as well, the PRIOR above is a better
% choice, as LOW_PRIOR might enforce too much the diffusion coefficient
% towards convergence in known values
LOW_PRIOR = squeeze(nanmean(nanmean((LSQNONNEG_AVG_PROFILE(:,:,round(size(data.img,3)/2),:))))); 
LOW_PRIOR(DRange > 500e-3) = 0;
LOW_PRIOR = LOW_PRIOR/sum(LOW_PRIOR);

% CHANGE
% disp('Different assignment')
if(parameters.use_low_prior == 1)
    PRIOR = LOW_PRIOR; 
end

% Heuristics to subdivide the number of classes, also taking into account
% restrictions on their own "_K"
[ClassAssignment,FinalMergedGaussians,FinalMergedGaussians_D,FinalMergedGaussians_f,OFinalMergedGaussians] = SFGMMFit(DRange,Dictionary,GMMDictionary,LOW_PRIOR,0.15,1);
[ClassAssignment_K,FinalMergedGaussians_K,FinalMergedGaussians_K_D,FinalMergedGaussians_K_f] = SFGMMFit(DRange,Dictionary,GMMDictionary,LOW_PRIOR,0.25,0);

LambdaGrid = parameters.LambdaGrid; % The regularization parameter (L2)
IX = 1;

for x=1:size(data.img,1)
disp(['Total: ' num2str(x/size(data.img,1)*100) '%']);
for y=1:size(data.img,2)
for z=1:size(data.img,3)
    if(data.mask(x,y,z) == 0 || isnan(data.mask(x,y,z)))
        continue
    end
    NoisySimulatedSignal = double(squeeze(data.img(x,y,z,:)));
    if(sum(NoisySimulatedSignal == 0) > 1)
        continue
    end
    
    NoisySimulatedSignal = NoisySimulatedSignal/mean(NoisySimulatedSignal(data.bvals <= min_b));
    
    goodpoints = squeeze(OUTLIERS(x,y,z,:))==1;
   
%     unique_bvals.bvals = unique(round(data.bvals(goodpoints)));
%     UniqueDictionary = DW_BuildIsotropicDictionary(unique_bvals,DRange);

    L2LSQNONNEG_DECONV = DW_RegularizedDeconv(Dictionary(goodpoints,:),NoisySimulatedSignal(goodpoints),op_lp,LambdaGrid(IX),PRIOR);
    L2LSQNONNEG_PRIOR_AVG_PROFILE(x,y,z,:) = L2LSQNONNEG_DECONV;
end
end
end

% save([prefix suffix],'-v7.3')
% Create fractional / diffusion maps and save

MAX_COMPONENTS = max(ClassAssignment);
IntegrationRanges = cell(max(ClassAssignment),1);
for ij=1:max(ClassAssignment)
   IntegrationRanges{ij} = find(ClassAssignment == ij);
end

MAX_COMPONENTS_K = max(ClassAssignment_K);
IntegrationRanges_K = cell(max(ClassAssignment_K),1);
for ij=1:max(ClassAssignment_K)
   IntegrationRanges_K{ij} = find(ClassAssignment_K == ij);
end

do_clean = 1;


TARGET = L2LSQNONNEG_PRIOR_AVG_PROFILE;

clean_data = DW_LoadData(data_name,bvec_name,bval_name,mask_name,0,'',1);
clean_data.bvals = round(clean_data.bvals);
GP = clean_data.bvals >= parameters.min_bval & clean_data.bvals < parameters.max_bval;
clean_data.img = clean_data.img(:,:,:,GP);
clean_data.bvals = clean_data.bvals(GP);
clean_data.bvecs = clean_data.bvecs(GP,:);
clean_data.hdr.dime.dim(5) = sum(GP);

NGnii = clean_data;
NGnii.img = NGnii.img(:,:,:,1);
NGnii.hdr.dime.dim(1) = 3;
NGnii.hdr.dime.dim(5) = 1;
NGnii.hdr.dime.bitpix = 32;
NGnii.hdr.dime.datatype = 16;

FractionsNii = clean_data;
FractionsNii.img = zeros([size(FractionsNii.img(:,:,:,1)) MAX_COMPONENTS]);
FractionsNii.hdr.dime.dim(5) = MAX_COMPONENTS;
FractionsNii.hdr.dime.bitpix = 32;
FractionsNii.hdr.dime.datatype = 16;

LSQNONNEG_CL_f = zeros([size(data.img(:,:,:,1)) MAX_COMPONENTS]);
LSQNONNEG_CL_D_mu = zeros([size(data.img(:,:,:,1)) MAX_COMPONENTS]);
LSQNONNEG_CL_D_sd = zeros([size(data.img(:,:,:,1)) MAX_COMPONENTS]);
NG = zeros([size(data.img(:,:,:,1)) MAX_COMPONENTS_K]);

for x=1:size(data.img,1)
disp(['Total: ' num2str(x/size(data.img,1)*100) '%']);
for y=1:size(data.img,2)
for z=1:size(data.img,3) 
    if(data.mask(x,y,z) == 0 || isnan(data.mask(x,y,z)))
        continue
    end
        deconv = squeeze(TARGET(x,y,z,:))';
        deconv(DRange > 500e-3) = 0;
        deconv(end) = 0;
        [LSQNONNEG_CL_D_mu(x,y,z,:),LSQNONNEG_CL_D_sd(x,y,z,:),LSQNONNEG_CL_f(x,y,z,:)] = SpectrumStatisticsInIntervals(DRange,deconv,IntegrationRanges);            
        [~,~,NG(x,y,z,:)] = SpectrumStatisticsInIntervals(DRange,deconv,IntegrationRanges_K);            
        NG(x,y,z,:) = NG(x,y,z,:)/sum(NG(x,y,z,:),4);
        LSQNONNEG_CL_f(x,y,z,:) = LSQNONNEG_CL_f(x,y,z,:)/sum(LSQNONNEG_CL_f(x,y,z,:),4);

end
end
end
if(do_clean == 1)
    NGnii.img = NG(:,:,:,1);
    save_untouch_nii(NGnii,[parameters.output_prefix '_NG.nii']);
    FractionsNii.img = LSQNONNEG_CL_f;
    save_untouch_nii(FractionsNii,[parameters.output_prefix '_f.nii']);
    FractionsNii.img = LSQNONNEG_CL_D_mu;
    save_untouch_nii(FractionsNii,[parameters.output_prefix '_D.nii']);
end

end

% Helper function for the spectral fit
function [Mean,Sd,f] = SpectrumStatisticsInIntervals(DRange,spectrum,Intervals)
    Mean = zeros(length(Intervals),1);
    Sd = zeros(length(Intervals),1);
    f = zeros(length(Intervals),1);
    
    for ij=1:length(Intervals)
       f(ij) = sum(spectrum(Intervals{ij}));
       Mean(ij) = sum(DRange(Intervals{ij}).*spectrum(Intervals{ij}))/sum(spectrum(Intervals{ij})); 
       M = length(find(spectrum(Intervals{ij}) > 0));
       if(M > 0)
           Sd(ij) = sqrt(sum(spectrum(Intervals{ij}).*(DRange(Intervals{ij})-Mean(ij)).^2)/((M-1)/M*sum(spectrum(Intervals{ij}))));
       end
    end
end

% Helper function for the spectral fit
function [musDictionary,mu_sd_couples] = DW_BuildGMMDictionaryWithCR(DMuRange,ratios)
NSignals = length(DMuRange);
musDictionary = zeros(NSignals,NSignals*length(ratios));
mu_sd_couples = zeros(length(DMuRange)*length(ratios),2);
index = 1;
for j=1:NSignals
    for k=1:length(ratios)
        musDictionary(:,(j-1)*length(ratios)+k) = Gaussian1D(DMuRange',DMuRange(j),DMuRange(j)/ratios(k));
        mu_sd_couples(index,:) = [DMuRange(j) DMuRange(j)/ratios(k)];
        index = index+1;
    end
end
end

% Helper function for the spectral fit
function Dictionary = DW_BuildIsotropicDictionary(data,DRange)

    Dictionary = zeros(length(data.bvals),length(DRange));

    for j=1:length(DRange)
        Dictionary(:,j) = exp(-data.bvals*DRange(j));
    end

end

% Helper function for the spectral fit
function data = DW_LoadData(data_name,bvec_name,bval_name,mask_name,load_nii_instead,noise_map,round_small_bs_to_zero,force_double)
if(exist('load_nii_instead','var') > 0 && load_nii_instead == 1)
    data = load_nii(data_name);
    disp('Loading with load_nii');
else
    data = load_untouch_nii(data_name);
    disp('Loading with load_untouch_nii');
end

data.bvecs = load(bvec_name);
data.bvals = load(bval_name);
data.bvecs = data.bvecs';
data.bvals = data.bvals';
data.img(data.img<0) = 0;
if(data.hdr.dime.scl_slope == 0)
    data.hdr.dime.scl_slope = 1;
end

if(exist('force_double','var') > 0 && force_double == 1)
    data.img = double(data.img)*data.hdr.dime.scl_slope+data.hdr.dime.scl_inter;
else
    data.img = single(data.img)*data.hdr.dime.scl_slope+data.hdr.dime.scl_inter;
end

if(exist('round_small_bs_to_zero','var') > 0 && round_small_bs_to_zero > 0)
    disp('Rounding small bs to zero');
    data.bvals(data.bvals < 1) = 0;
end

if(exist('mask_name','var') > 0 && ~isempty(mask_name))
    mask = load_untouch_nii(mask_name);
    data.mask = mask.img;
else
    data.mask = data.img(:,:,:,find(data.bvals<1,1))>0;
end
if(size(data.img,4) > length(data.bvals))
%     WE HAVE NOISE MAP AT THE END
    disp('Noisemap found');
    data.noisemap = data.img(:,:,:,end);
    data.img = data.img(:,:,:,1:end-1);
    data.hdr.dime.dim(5) = size(data.img,4);
else
    if(exist('noise_map', 'var') > 0 && ~isempty(noise_map))
        disp('Loading noise map');
        try
        if(exist('load_nii_instead','var') > 0 && load_nii_instead == 1)
            nmap = load_nii(noise_map);
        else
            nmap = load_untouch_nii(noise_map);
        end
        data.noisemap = nmap.img;
        catch
            disp('Error loading noisemap, please check ');
        end
    else
        %INDEED CHECK IF A NOISE MAP EXIST
        base_name = strsplit(data_name,'.');
        if(exist([base_name{1} '_noisemap.nii.gz'],'var') > 0)
           disp('Found a noisemap matching dataname, loading'); 
            if(exist('load_nii_instead','var') > 0 && load_nii_instead == 1)
                nmap = load_nii([base_name{1} '_noisemap.nii.gz']);
            else
                nmap = load_untouch_nii([base_name{1} '_noisemap.nii.gz']);
            end
            data.noisemap = nmap.img;
        end
    end
end
end

% Helper function for the spectral fit
function DW_SaveData(data,out_name,niigz)
    if(nargin < 3)
        niigz = 1;
    end
    if(niigz == 1)
        ext = '.nii.gz';
    else
        ext = '.nii';
    end
    data.img = data.img/data.hdr.dime.scl_slope-data.hdr.dime.scl_inter;
    save_untouch_nii(data,[out_name ext]);
    if(isfield(data,'bvals'))
        bvals = data.bvals';
        bvecs = data.bvecs';
        save([out_name '.bvec'],'bvecs','-ascii');
        save([out_name '.bval'],'bvals','-ascii');
    end
    if(isfield(data,'mask') > 0)
       mask.hdr = data.hdr;
       mask.untouch = 1;
       mask.hdr.dime.dim(1) = 3;
       mask.hdr.dime.dim(5) = 1;
       mask.img = zeros(mask.hdr.dime.dim(2:4));
       mask.img = data.mask;
       save_untouch_nii(mask,[out_name '_mask' ext]);
    end
    if(isfield(data,'noisemap') > 0)
       mask.hdr = data.hdr;
       mask.untouch = 1;
       mask.hdr.dime.dim(1) = 3;
       mask.hdr.dime.dim(5) = 1;
       mask.img = zeros(mask.hdr.dime.dim(2:4));
       mask.img = data.noisemap;
       save_untouch_nii(mask,[out_name '_noisemap' ext]);
    end
end

% Helper function for the spectral fit
function [Deconv_1,Deconv_1_clean] = DW_RegularizedDeconv(Dictionary_1,OneBigVoxelFull_1,options,lambda,x0)
if(nargin < 3)
    options = optimset('TolX',1e-12,'TolFun',1e-12);
end

% base_deconv = lsqnonneg(Dictionary_1,OneBigVoxelFull_1,options);
% Amp = sum(base_deconv);
Amp = 1;

C = Amp*lambda*eye(size(Dictionary_1,2),size(Dictionary_1,2));
C = C(1:size(Dictionary_1,2),1:size(Dictionary_1,2));
Dictionary_1_padded = [Dictionary_1; C];
OneBigVoxelFull_1_tmp = zeros(size(Dictionary_1_padded,1),1);
OneBigVoxelFull_1_tmp(1:length(OneBigVoxelFull_1)) = OneBigVoxelFull_1;
if(exist('x0','var') > 0)
   OneBigVoxelFull_1_tmp(length(OneBigVoxelFull_1)+1:end) = x0; 
end

Deconv_1 = lsqnonneg(Dictionary_1_padded,double(OneBigVoxelFull_1_tmp),options);

if(nargout > 1)
   Deconv_1_clean = zeros(size(Deconv_1));
   idx = Deconv_1>0;
   Deconv_1_clean(idx) = lsqnonneg(Dictionary_1(:,idx),OneBigVoxelFull_1,options);
end

end

% Helper function for the spectral fit
function [LSQNONNEG_CL_DECONV,goodpoints] = DW_RobustDeconvFinal(Dictionary,NoisySimulatedSignal,kappa,options)

goodpoints = ones(size(NoisySimulatedSignal)) == 1;%squeeze(OUTLIERS(x,y,z,:))==0;
W = diag(NoisySimulatedSignal);
done = 0;
iters = 0;
if(nargin < 3)
    kappa = 3;
    options = optimset('TolX',1e-2);
end

while done == 0
    iters = iters+1;
    new_outliers_hes = 0;
    new_outliers_hos = 0;
    for iter=1:4
        LSQNONNEG_CL_DECONV = lsqnonneg(Dictionary(goodpoints,:),NoisySimulatedSignal(goodpoints),options);
        residuals = NoisySimulatedSignal-Dictionary*LSQNONNEG_CL_DECONV;
        sd = mad(residuals)*1.4826;
        bad_guys = abs(residuals-mean(residuals))>kappa*sd;
        if(sum(bad_guys(goodpoints == 1)==1) > 0)
            new_outliers_hes = 1;
        end
        goodpoints = goodpoints & bad_guys == 0;
        %             disp(['Outliers (HeS): ' num2str(length(find(goodpoints==0)))]);
    end
    for iter=1:4
        LSQNONNEG_CL_DECONV = lsqnonneg(W(goodpoints,goodpoints)*Dictionary(goodpoints,:),W(goodpoints,goodpoints)*NoisySimulatedSignal(goodpoints),options);
        W = diag((Dictionary*LSQNONNEG_CL_DECONV));
        residuals = NoisySimulatedSignal-Dictionary*LSQNONNEG_CL_DECONV;
        sd = mad(residuals)*1.4826;
        bad_guys = abs(residuals-mean(residuals))>kappa*sd;
        if(sum(bad_guys(goodpoints == 1)==1) > 0)
            new_outliers_hos = 1;
        end
        goodpoints = goodpoints & bad_guys == 0;
        %             disp(['Outliers (HoS): ' num2str(length(find(goodpoints==0)))]);
    end
    if(new_outliers_hos == 0 && new_outliers_hes == 0)
        done = 1;
    end
    if(iters >= 100)
        done = 1;
    end
end

LSQNONNEG_CL_DECONV = lsqnonneg(Dictionary(goodpoints,:),NoisySimulatedSignal(goodpoints),options);

end

% Helper function for the spectral fit
function out = Gaussian1D(x,mu,sigma)
%     out = (1/(sqrt(2*pi)*sigma))*exp(-0.5*(mu-x).^2/sigma^2);
    out = exp(-0.5*(mu-x).^2/sigma^2);
end

% Helper function for the spectral fit
function [ClassAssignment,FinalMergedGaussians,FinalMergedGaussians_D,FinalMergedGaussians_f,OFinalMergedGaussians] = SFGMMFit(DRange,Dictionary,GMMDictionary,U2,MinOverlap2Merge,FuseComp12)
op_e8 = optimset('TolX',1e-8);
if(exist('MinOverlap2Merge','var') < 1)
    MinOverlap2Merge = 0.15;
end
if(exist('FuseComp12','var') < 1)
    FuseComp12 = 1;
end

REC = Dictionary*U2;%(U/sum(U));
IDEC = lsqnonneg(Dictionary*GMMDictionary,REC,op_e8);
CP = find(IDEC);
IndividualGaussians = zeros(length(DRange),length(CP));
for ij=1:length(CP)
   DCopy = zeros(size(IDEC));
   DCopy(CP(ij)) = IDEC(CP(ij));
   IndividualGaussians(:,ij) = GMMDictionary*DCopy; 
   
   Area = sum(IndividualGaussians(:,ij));
   Ak = cumsum(IndividualGaussians(:,ij));
   IndividualGaussians(1:find(Ak >= 0.025*Area,1),ij) = 0;
   IndividualGaussians(find(Ak >= 0.975*Area,1):end,ij) = 0;

end
Overlap = zeros(length(CP),length(CP));
for ij=1:length(CP)
    for ik=1:length(CP)
%         Overlap(ij,ik) = sum(min(IndividualGaussians(:,ij),IndividualGaussians(:,ik)))/sum(max(IndividualGaussians(:,ij),IndividualGaussians(:,ik)));
        Overlap(ij,ik) = sum(min(IndividualGaussians(:,ij),IndividualGaussians(:,ik)))/min(sum(IndividualGaussians(:,ij)),sum(IndividualGaussians(:,ik)));
    end
end
MergeComponents = 1:length(CP);
no_change = 1;
    nof_changes = 0;
    for ij=1:length(CP)
        for ik=1:length(CP)
            if(MergeComponents(ik) ~= ij && Overlap(ij,ik) > MinOverlap2Merge && Overlap(ij,ik) >= max(Overlap(ik,[1:ik-1 ik+1:end])))
               if(MergeComponents(ij) == ij)
                   MergeComponents(ik) = ij;
               else
                   MergeComponents(ik) = MergeComponents(ij);
               end
               no_change = 1;
               nof_changes = nof_changes+1;
            end
        end
    end
    disp(nof_changes)
    if(nof_changes == 0)
        no_change = 0;
    end
    
UQ = unique(MergeComponents);
MergedGaussians = zeros(length(DRange),length(UQ));
MergedGaussians_D = zeros(length(UQ),1);
MergedGaussians_f = zeros(length(UQ),1);
for ij=1:length(UQ)
   MergedGaussians(:,ij) = sum(IndividualGaussians(:,MergeComponents == UQ(ij)),2); 
   [MergedGaussians_D(ij),~,MergedGaussians_f(ij)] = SpectrumStatisticsInIntervals(DRange,MergedGaussians(:,ij)',{1:length(DRange)});
end

% Simplification rules
FinalMergeComponents = 1:length(UQ);
if(FuseComp12 == 1)
    FinalMergeComponents(find(MergedGaussians_D < 0.3e-3,1)) = find(MergedGaussians_D > 0.3e-3,1);
end
FinalMergeComponents(MergedGaussians_D > 6e-3) = find(MergedGaussians_D > 6e-3,1);

UQ = unique(FinalMergeComponents);
FinalMergedGaussians = zeros(length(DRange),length(UQ));
OFinalMergedGaussians = FinalMergedGaussians;
FinalMergedGaussians_D = zeros(length(UQ),1);
FinalMergedGaussians_f = zeros(length(UQ),1);
for ij=1:length(UQ)
   FinalMergedGaussians(:,ij) = sum(MergedGaussians(:,FinalMergeComponents == UQ(ij)),2); 
   OFinalMergedGaussians(:,ij) = FinalMergedGaussians(:,ij);
   Area = sum(FinalMergedGaussians(:,ij));
   Ak = cumsum(FinalMergedGaussians(:,ij));
   FinalMergedGaussians(1:find(Ak >= 0.025*Area,1),ij) = 0;
   FinalMergedGaussians(find(Ak >= 0.975*Area,1):end,ij) = 0;
%    FinalMergedGaussians(FinalMergedGaussians(:,ij) < 0.02*max(FinalMergedGaussians(:,ij)),ij) = 0;
   [FinalMergedGaussians_D(ij),~,FinalMergedGaussians_f(ij)] = SpectrumStatisticsInIntervals(DRange,FinalMergedGaussians(:,ij)',{1:length(DRange)});
end

ClassAssignment = zeros(length(DRange),1);

DIX = zeros(length(FinalMergedGaussians_D),1);
for ij=1:length(DIX)
   DIX(ij) = find(DRange >= FinalMergedGaussians_D(ij),1); 
end
for ij=1:length(ClassAssignment)
   [m,IX] = max(FinalMergedGaussians(ij,:));
   if(m ~= 0)
       rM = max(ClassAssignment(1:ij-1));
       if(~isempty(rM))
           IX = max(IX,rM);
       end
       ClassAssignment(ij) = IX;
   else
%        [~,IX] = min(abs(DRange(ij)-FinalMergedGaussians_D));
       [~,IX] = min(abs(ij-DIX));
       rM = max(ClassAssignment(1:ij-1));
       if(~isempty(rM))
           IX = max(IX(1),rM);
       end
       ClassAssignment(ij) = IX(1);
   end
end
end

% Geometrical average per shell/unique b-val
function s1 = QuickGeoAverage(data)
   data.bvals = round(data.bvals);
   ub1 = unique(data.bvals);

   s1 = data;
   s1.img = zeros([size(data.img(:,:,:,1)),length(ub1)]);

   s1.bvals = ub1;
   s1.bvecs = [];

   for j=1:length(ub1)
      vol = single(data.img(:,:,:,data.bvals==ub1(j)));
      s1.img(:,:,:,j) = geomean(vol,4);
   end
        
end

% Segmented IVIM + tissue (isotropic) fit
function [Dhigh,Dlow,f] = SegmentedIVIMFit2Comp(data,parameters)
    data = QuickGeoAverage(data);
    IX_high = data.bvals >= parameters.bthresh & data.bvals <= parameters.max_bval;
    IX_low = data.bvals >= parameters.min_bval & data.bvals < parameters.bthresh;

    siz = size(data.img);
    data.img = log(reshape(data.img,siz(1)*siz(2)*siz(3),siz(4)));
    Xh = [-data.bvals(IX_high) ones(sum(IX_high==1),1)];
    HighD = Xh\data.img(:,IX_high)';

    Dhigh.img = reshape(HighD(1,:),siz(1:3));
    Dhigh.VD = data.VD;
    Dhigh.img(data.mask == 0) = 0;

    Xl =[-data.bvals(IX_low) ones(sum(IX_low==1),1)];
    LowD = Xl\data.img(:,IX_low)';
    Dlow.img = reshape(LowD(1,:),siz(1:3));
    Dlow.VD = data.VD;
    Dlow.img(data.mask == 0) = 0;

    data.img = exp(data.img);
    f.img = abs(data.img(:,1) - exp(HighD(2,:)'))./data.img(:,1);
    f.img = reshape(f.img,siz(1:3));
    f.VD = data.VD;
    f.img(data.mask == 0) = 0;
end

% FW + D model
function [out,S] = FW2CompartmentsModel(x0,params)
    S = x0(1)*((1-x0(2))*exp(-params.b*x0(3))+x0(2)*exp(-params.b*3e-3));
    out = params.S - S;
end

% FW + Tissue (isotropic) non-linear fit
function [Dhigh,f] = FWFit2Comp(data,parameters)
    data = QuickGeoAverage(data);
    IX_high = data.bvals >= parameters.min_bval & data.bvals <= parameters.max_bval;

    x0 = [1 0.05 0.7e-3];
    lb = [0 0 0];
    ub = [1 1 2e-3];
    
    Dhigh.img = zeros(size(data.img(:,:,:,1)));
    Dhigh.VD = data.VD;
    f.img = zeros(size(data.img(:,:,:,1)));
    f.VD = data.VD;
    
    siz = size(data.img);
    
    data.img = (reshape(data.img,siz(1)*siz(2)*siz(3),siz(4)));
    data.img = permute(data.img,[2 1]);
    points2fit = find(data.mask > 0);

    options = optimset('TolX',1e-3,'TolFun',1e-3,'Display','off');
    params.b = data.bvals(IX_high);
    tic
    for pi=1:length(points2fit)
        if(mod(pi,1000) == 0)
            disp(['Free-water fit: ' num2str(pi/length(points2fit)*100) '%']);
        end
        params.S = data.img(IX_high,points2fit(pi));
        x0(1) = params.S(1);
        ub(1) = 2*params.S(1);
        p = lsqnonlin(@FW2CompartmentsModel,x0,lb,ub,options,params);
        Dhigh.img(points2fit(pi)) = p(3);
        f.img(points2fit(pi)) = p(2);
    end
    toc
end

% IVIM + FW + D model
function [out,S] = IVIMFW3CompartmentsModel(x0,params)
    S0 = x0(1);
    fivim = x0(2);
    ffw = x0(3);
    D = x0(4);
    DStar = x0(5);
    S = S0*((1-fivim)*((1-ffw)*exp(-params.b*D)+ffw*exp(-params.b*3e-3)) +...
        fivim*exp(-params.b*DStar));
    out = params.S - S;
end

% IVIM + FW + Tissue (isotropic) non-linear fit
function [Dhigh,DStar,ffw,fivim] = IVIMFWFit3Comp(data,parameters)
    data = QuickGeoAverage(data);
    IX_high = data.bvals >= parameters.min_bval & data.bvals <= parameters.max_bval;

    x0 = [1 0.05 0.05 0.7e-3 10e-3];
    lb = [0 0 0 0 6e-3];
    ub = [1 1 1 2e-3 500e-3];
    
    Dhigh.img = zeros(size(data.img(:,:,:,1)));
    Dhigh.VD = data.VD;
    ffw.img = zeros(size(data.img(:,:,:,1)));
    ffw.VD = data.VD;
    fivim.img = zeros(size(data.img(:,:,:,1)));
    fivim.VD = data.VD;
    DStar.img = zeros(size(data.img(:,:,:,1)));
    DStar.VD = data.VD;
    
    siz = size(data.img);
    
    data.img = (reshape(data.img,siz(1)*siz(2)*siz(3),siz(4)));
    data.img = permute(data.img,[2 1]);
    points2fit = find(data.mask > 0);

    options = optimset('TolX',1e-3,'TolFun',1e-3,'Display','off');
    params.b = data.bvals(IX_high);
    tic
    for pi=1:length(points2fit)
        if(mod(pi,1000) == 0)
            disp(['IVIM + free-water fit: ' num2str(pi/length(points2fit)*100) '%']);
        end
        params.S = data.img(IX_high,points2fit(pi));
        x0(1) = params.S(1);
        ub(1) = 2*params.S(1);
        p = lsqnonlin(@IVIMFW3CompartmentsModel,x0,lb,ub,options,params);
        Dhigh.img(points2fit(pi)) = p(4);
        fivim.img(points2fit(pi)) = p(2);
        ffw.img(points2fit(pi)) = p(3);
        DStar.img(points2fit(pi)) = p(5);
    end
    toc
end

