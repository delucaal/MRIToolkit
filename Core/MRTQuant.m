

%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

% This class implements methods for the MRIToolkit diffusion pipeline
% Differently from other classes (as EDTI), the input / output of this
% class are Niftis and other exchangeble formats (in addition to
% ExploreDTI-compatible .mat files)
classdef MRTQuant < handle
    methods(Static)
        
        function data = LoadNifti(data_file,apply_intensity_scale)
        % Loads a .nii(.gz) and returns it into a data structure with
        % fields: img (the data matrix) and VD (voxel dimension)
        % input arguments:
        % data_file: the .nii file to load
        % apply_intensity_scale: use the slope/intercept specified in the
        % nifti header
            if(nargin < 2)
                apply_intensity_scale = 0;
            end
            [I,VD,~,hdr] = EDTI_Library.E_DTI_read_nifti_file(data_file);
            data.img = I;
            data.VD = VD;
            data.hdr = hdr;
            
            if(apply_intensity_scale == 1 && hdr.dime.scl_slope ~= 0 && ...
                    isfinite(hdr.dime.scl_slope) && isfinite(hdr.dime.scl_inter))
                data.img = single(data.img)*hdr.dime.scl_slope+hdr.dime.scl_inter;
            end
        end
        
        function WriteNifti(data, file_name)
        % Writes a data matrix to nifti as specified in file_name. data
        % must be a struct with fields img (the matrix) and VD (voxel
        % dimension)
            EDTI_Library.E_DTI_write_nifti_file(data.img, data.VD, file_name);
        end
        
        function ConformSpatialDimensions(varargin)
        % Ensures consistency with the coordinate systems of MRIToolkit and
        % ExploreDTI. Input arguments:
        % nii_file: the input file
        % output: the output file
            if( isempty(varargin))
                my_help('MRTQuant.ConformSpatialDimensions');
                return;
            end

            coptions = varargin;
            nii_file = GiveValueForName(coptions,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory parameter nii_file');
            end
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                error('Missing mandatory parameter output');
            end
            
            data_in = nii_file;
            data_out = output;

            hdr = load_untouch_nii(data_in);
            smat = [hdr.hdr.hist.srow_x(1:3);
                    hdr.hdr.hist.srow_y(1:3);
                    hdr.hdr.hist.srow_z(1:3)];
            permute_order = [1 2 3];
            sign_order = [1 1 1];
            for ij=1:size(smat,2)
                [~,IX] = max(abs(smat(:,ij)));
                permute_order(ij) = IX;
                if(smat(IX,ij) < 0)
                    sign_order(IX) = -1;
                end
            end

            out_data.img = single(hdr.img);
            out_data.img = permute(out_data.img,[permute_order 4]);
            for ij=1:length(sign_order)
               if(sign_order(ij) == -1)
                   out_data.img = flip(out_data.img,ij);
               end
            end

            out_data.VD = hdr.hdr.dime.pixdim(2:4);
            out_data.VD = out_data.VD(permute_order(1:3));
            out_data.VD = out_data.VD([2 1 3]);
            out_data.img = permute(out_data.img,[2 1 3 4]);
            out_data.img = flip(out_data.img,1);
            out_data.img = flip(out_data.img,2);

            MRTQuant.WriteNifti(out_data,data_out);

        end
        
        function SetB0Val(newval)
        %  Set the b-value threshold in s/mm2 to be considered a b=0s/mm2.
        % This is useful for instance when working with HCP data (set it to
        % 70).
            global MRIToolkit;
            MRIToolkit.min_bval_as_b0 = newval;
        end
        
        function bmat = b_Matrix_from_bval_bvec(varargin)
        % Computes bmat (b-matrix) from bvals and bvecs and writes it
        % to a .txt file with the same name. Specify only the
        % bval file (or variable) as input. The bvec must then be in the same folder with
        % the same naming. txt_file is the destination file of the b-matrix
        % (optional, otherwise the file is created in the same directory)
        % Arguments:
        % bval_file: the input .bval file
        % bvec_file: the input .bvec file - optional. If not specified a
        % file with the same name of bval_file and extension .bvec will be
        % searched.
        % or, alternatively
        % bval: the input bvals as variable
        % bvec: the input bvecs as variable
        
        % output: the output b-matrix .txt file
            if( isempty(varargin))
                my_help('MRTQuant.b_Matrix_from_bval_bvec');
                bmat = [];
                return;
            end
            coptions = varargin;
            bval = GiveValueForName(coptions,'bval');
            bvec = GiveValueForName(coptions,'bvec');
            if(isempty(bval) || isempty(bvec))
                bval_file = GiveValueForName(coptions,'bval_file');
                if(isempty(bval_file))
                    error('Missing mandatory parameter bval_file');
                end
                bvec_file = GiveValueForName(coptions,'bvec_file');
                if(isempty(bvec_file))
                   bvec_file = [bval_file(1:end-4) 'bvec'];
                end
                txt_file = GiveValueForName(coptions,'output');
                if(isempty(txt_file))
                    error('Missing mandatory parameter output');
                end
                EDTI_Library.E_DTI_convert_nii_dic_2_txt_exe2(bval_file,bvec_file,txt_file);
            else
                    g = bvec;
                    bmat = repmat(bval,[1 6]).*[g(:,1).^2 2*g(:,1).*g(:,2) 2*g(:,1).*g(:,3) ...
                        g(:,2).^2 2*g(:,2).*g(:,3) g(:,3).^2];
                    return
            end
            if(nargout > 0)
                if(isempty(txt_file))
                    bmat = load([bval_file(1:end-5) '.txt']);
                else
                    bmat = load(txt_file);
                end
            end
        end
        
        function [bvals,grad] = bval_bvec_from_b_Matrix(b)
        % Computes bvals and bvecs for a given b (b-matrix)
            bvals = sum(b(:,[1 4 6]),2);
            
            BSign = sign(sign(b(:,1:3)) + 0.0001);   % Adding 0.0001 avoids getting zeros here
            
            grad = BSign .* sqrt(b(:,[1 4 6]));
            grad_n = sqrt(sum(grad.*grad,2));
            grad = grad./[grad_n grad_n grad_n];
            grad(isnan(grad)) = 0;
        end
        
        function mask = EnhancedSpatialMasking(DWI, tuning_param, morphological_size)
        % Computes a mask for dMRI data using combinations of thresholding
        % and image-processing (erosions)
            mask = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced_IND(DWI, tuning_param, morphological_size);
        end
        
        function DataMat2Nii_bval_bvec(file_in,prefix_out)
        % Exports the data (.nii) + bvals and bvecs for a given file_in
        % (must be an ExploreDTI-like .mat file). The export will use the
        % naming specified in prefix_out
            MAT = load(file_in,'DWI','b','VDims');
            [bval, grad] =  EDTI_Library.E_DTI_GetGradientsandBval_SC(MAT.b, 0);
            MAT.DWI = EDTI_Library.E_DTI_DWI_cell2mat(MAT.DWI);
            
            EDTI_Library.E_DTI_write_nifti_file(MAT.DWI,MAT.VDims,[prefix_out '.nii']);
            
            bval_file = fopen([prefix_out '.bval'],'wt');
            for ij=1:length(bval)
                fprintf(bval_file,'%f ',bval(ij));
            end
            fclose(bval_file);
            
            bvec_file = fopen([prefix_out '.bvec'],'wt');
            for ij=1:size(grad,2)
                for ik=1:size(grad,1)
                    fprintf(bvec_file,'%f ',grad(ik,ij));
                end
                fprintf(bvec_file,'%s',newline);
            end
            fclose(bvec_file);
        end
        
        function PerformSignalDriftCorrection(varargin)
        % Interface to the signal drift correction. Arguments:
        % nii_file: the input .nii file
        % output: the output .nii file
        % target_bval: the target b-value used for the correction, default is 0
        % target_bval_tol: tolerance around the target b-value - default is 1
        % pol_degree: Degree of the polynomial used for correction (1-3), default is 2
        % masking: try automatic masking on the target volume. default is 0
            if(isempty(varargin))
                my_help('MRTQuant.PerformSignalDriftCorrection');
                return;
            end
            
            json.CallFunction = 'MRTQuant.PerformSignalDriftCorrection';
            json.Description = my_help('MRTQuant.PerformSignalDriftCorrection');

            coptions = varargin;
            par.bvalC = 0;
            tgt_bval = GiveValueForName(coptions,'target_bval');
            if(~isempty(tgt_bval))
                par.bvalC = (tgt_bval);
            end
            json.target_bval = tgt_bval;
            
            par.bv_thresh=1;
            tgt_thresh = GiveValueForName(coptions,'target_bval_tol');
            if(~isempty(tgt_thresh))
                par.bv_thresh = (tgt_thresh);
            end
            json.bv_thresh = tgt_thresh;
            
            par.method = 2;
            tgt_method = GiveValueForName(coptions,'pol_degree');
            if(~isempty(tgt_method))
                par.method = (tgt_method);
            end
            json.method = tgt_method;
            
            par.masking.do_it = 0;
            tgt_masking = GiveValueForName(coptions,'masking');
            if(~isempty(tgt_masking))
                par.tgt_masking = (tgt_masking);
            end
            json.masking = tgt_masking;
            
            par.show_summ_plot = 1;
            par.suff = '_sdc';
            par.masking.p1 = 5;
            par.masking.p2 = 1;
            
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Missing mandatory argument nii_file');
            end 
            
            file_out = GiveValueForName(coptions,'output');
            if(isempty(file_out))
                error('Missing mandatory argument output');
            end 
            
            param{1} = par;
            param{1}.f_in_nii = file_in;
            param{1}.f_in_txt = [file_in(1:end-4) '.txt'];
            param{1}.f_out_nii = file_out;
            param{1}.f_out_txt = [file_out(1:end-4) '.txt'];
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';
            
            EDTI_Library.E_DTI_signal_drift_correction(param{1});
            
            NiftiIO_basic.WriteJSONDescription('output',file_out(1:end-4),'props',json);
        end
        
        function FlipPermuteSpatialDimensions(varargin)
        % Interface to flip / permuting spatial dimensions (Ensures
        % compliance of the .nii files with the ExploreDTI convention)
        % Accordingly, the header is discarded and only the rotation matrix
        % retained. Arguments:
        % nii_file: The file to flip permute.
        % output: The flip/permuted file.
        % permute: default is [1 2 3]
        % flip: default is [0 0 0]. Setting 1 will flip the corresponding
        %       dimension
            if(isempty(varargin))
                my_help('MRTQuant.FlipPermuteSpatialDimensions');
                return;
            end
            
            json.CallFunction = 'MRTQuant.FlipPermuteSpatialDimensions';
            json.Description = my_help('MRTQuant.FlipPermuteSpatialDimensions');
            
            coptions = varargin;
            params.suff = '';
            tgt_permute = GiveValueForName(coptions,'permute');
            tgt_flip = GiveValueForName(coptions,'flip');
            if(isempty(tgt_permute))
                params.permute = [1 2 3];
            else
                params.permute = tgt_permute;
            end
            if(isempty(tgt_flip))
                params.flip = [0 0 0];
            else
                params.flip = tgt_flip;
            end
            
            json.permute = params.permute;
            json.flip = params.permute;
            
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Missing mandatory parameter nii_file');
            end
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';

            file_out = GiveValueForName(coptions,'output');
            if(isempty(file_out))
                file_out = [file_in(1:end-4) '_FP.nii'];
            end
                        
            params.force_voxel_size = [];
            EDTI_Library.E_DTI_flip_permute_nii_file_exe(file_in,params,file_out);
                        
            if(exist([file_in(1:end-4) '.txt'],'file'))
                copyfile([file_in(1:end-4) '.txt'],[file_out(1:end-4) '.txt']);
            end
            
            NiftiIO_basic.WriteJSONDescription('output',file_out(1:end-4),'props',json);
        end

        function ResampleDataSpatially(varargin)
        % Inteface to resample 3D/4D volumes to a given resolution (Ensures
        % compliance of the .nii files with the ExploreDTI convention)
        % Accordingly, the header is discarded and only the rotation matrix
        % retained. Arguments:
        % nii_file: The file to flip permute.
        % output: The flip/permuted file.
        % res: the desired resolutio, e.g. [1 1 1]
            if(isempty(varargin))
                my_help('MRTQuant.ResampleDataSpatially');
                return;
            end

            json.CallFunction = 'MRTQuant.ResampleDataSpatially';
            json.Description = my_help('MRTQuant.ResampleDataSpatially');

            coptions = varargin;
            params.suff = '';
            tgt_res = GiveValueForName(coptions,'res');
            if(isempty(tgt_res))
                error('Missing mandatory parameter tgt_res');
            end            
            json.tgt_res = tgt_res;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Missing mandatory parameter nii_file');
            end
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';

            file_out = GiveValueForName(coptions,'output');
            if(isempty(file_out))
                file_out = [file_in(1:end-4) '_resampled.nii'];
            end
                        
            data = MRTQuant.LoadNifti(file_in);
            
            data.img = EDTI_Library.E_DTI_resample_nii_file(data.img, tgt_res, data.VD);
            data.VD = tgt_res;
            
            MRTQuant.WriteNifti(data,file_out);
            
            if(exist([file_in(1:end-4) '.txt'],'file'))
                copyfile([file_in(1:end-4) '.txt'],[file_out(1:end-4) '.txt']);
            end
            
            NiftiIO_basic.WriteJSONDescription('output',file_out(1:end-4),'props',json);
        end
        
        function PerformDTI_DKIFit(varargin)
        % Interface to perform the DTI/DKI fit. This will create a .mat
        % file in the ExploreDTI format with the same name of the provided
        % .nii file. Arguments:
        % nii_file: the diffusion MRI 4D .nii data
        % txt_file: The b-matrix .txt file. If empty, will look for a file
        %           named as the nii_file but with .txt extension.
        % dki: either 0 (default) or 1, enables the DKI fit
        % grad_perm: gradient permutations to ensure consistency with the
        %            spatial dimensions. 1[x y z] 2[y x z] 3[z y x]
        %            4[x z y] 5[y z x] 6[z x y]
        % grad_flip: sign flipping of the gradients to ensure consistency
        %            1[x y z] 2[-x y z] 3[x -y z] 4[x y -z]
        % fit_mode: ols, nls, wls, rekindle
        % dki_constraints: [min max] for the DKI diagonal values
        % mask_file: a .nii for masking
        % output: the output .mat file. If not specified a file with the
        % same name of the .nii but extension .mat is created
        % mk_curve: 1 or 0. Repeat the DKI fit with the "M-K curve technique"
            if(isempty(varargin))
                my_help('MRTQuant.PerformDTI_DKIFit');
                return;
            end                    
            global MRIToolkit;
            
            json.CallFunction = 'MRTQuant.PerformDTI_DKIFit';
            json.Description = my_help('MRTQuant.PerformDTI_DKIFit');
            
            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Need to specify the target file');
            end
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Quantification';
            
            txtfile = GiveValueForName(coptions,'txt_file');
            if(isempty(txtfile))
                txtfile = [file_in(1:end-4) '.txt'];
            end
            
            json.txt_file = txtfile;
            
            mk_curve_fit = GiveValueForName(coptions,'mk_curve');
            if(isempty(mk_curve_fit))
                mk_curve_fit = 0;
            end
            
            json.mk_curve_fit = mk_curve_fit;
            
            gmat = load(txtfile);
            bvals = sum(gmat(:,[1 4 6]),2);
            
            option_value = GiveValueForName(coptions,'dki');
            if(~isempty(option_value))
                if((option_value) > 0)
                    dki_fit = 1;
                else
                    dki_fit = 0;
                end
            else
                dki_fit = 0;
            end
            
            json.dki_fit = dki_fit;
            
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
                        
            option_value = GiveValueForName(coptions,'fit_mode');
            if(~isempty(option_value))
                fit_mode = option_value;
            else
                fit_mode = 'ols';
            end
            
            json.fit_mode = fit_mode;
            
            option_value = GiveValueForName(coptions,'dki_constraints');
            if(~isempty(option_value))
                V2 = strsplit(option_value,' ');
                if(length(V2) ~= 2)
                    error('Wrong format. dki_constraints=[Min Max]');
                end
                D1 = str2double(V2{1}(2:end));
                D2 = str2double(V2{2}(1:end-1));
                dki_constraints = [D1 D2];
            else
                dki_constraints = [-Inf Inf];
            end
            
            json.dki_constraints = dki_constraints;
            
            option_value = GiveValueForName(coptions,'output');
            if(~isempty(option_value))
                out_name = option_value;
            else
                out_name = [file_in(1:end-4) '.mat'];
            end
            
            option_value = GiveValueForName(coptions,'mask_file');
            if(~isempty(option_value))
                Mask_par = option_value;
                json.mask_file = Mask_par;
            else
                if(isfield(MRIToolkit,'EDTI_settings') && ...
                        isfield(MRIToolkit.EDTI_settings,'Mask_par'))
                    Mask_par.tune_NDWI = MRIToolkit.EDTI_settings.Mask_par.tune_NDWI;
                    Mask_par.tune_DWI = MRIToolkit.EDTI_settings.Mask_par.tune_DWI;
                    Mask_par.mfs = 7;
                else
                    Mask_par.tune_NDWI = 0.5;
                    Mask_par.tune_DWI = 0.5;
                    Mask_par.mfs = 7;
                end
            end
            
            if(isempty(perm) || isempty(flip))
               disp('Trying to automatically determine the coordinate systems. Warning: this works only with brain data!');
               while(true)
                  temp_file = fullfile(tempdir,['MRT_EDTI_' num2str(randi(500000)) '.mat']);
                  if(exist(temp_file,'file') < 1)
                      break
                  end
               end
            	EDTI_Library.E_DTI_model_fit(file_in,txtfile,temp_file,Mask_par,1,1, 'ols',0,[]);
               
               [flips, perms, correct, consistent] = EDTI_Library.E_DTI_Check_for_flip_perm_grads(temp_file);
               delete(temp_file)
               perm_list = [1 2 3;
                            2 1 3;
                            3 2 1;
                            1 3 2;
                            2 3 1;
                            3 1 2];
                sign_list = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1];
               if(isempty(flips))
                   if(~isdeployed)
                       warning('Failed to determine the coordinate systems. Please check your results / specify it manually.');
                   end
                   perm = 2;
                   flip = 2;
               else
                   for p=1:size(perm_list,1)
                      if(all(perm_list(p,:) == perms))
                          break
                      end
                   end
                   for f=1:size(sign_list,1)
                      if(all(sign_list(f,:) == flips))
                          break
                      end
                   end
                   perm = p;
                   flip = f;
               end
            end
            
            json.grad_flip = flip;
            json.grad_perm = perm;
            
            EDTI_Library.E_DTI_model_fit(file_in,txtfile,out_name,Mask_par,perm,flip, fit_mode,dki_fit, dki_constraints);
            
            if(mk_curve_fit == 1)
               disp('Performing M-K curve fit');
               mkcurve_fit_mat(out_name,out_name); 
            end
            
            NiftiIO_basic.WriteJSONDescription('output',out_name(1:end-4),'props',json);
        end
        
        function RefitDKIWithMKCurve(varargin)
           % This function re-fits an existing .mat file with the MK-curve
           % method (Zhang et al. 2019). Input arguments:
           % mat_file: The ExpoloreDTI-like .mat file
           % output: The output tracts (must be .mat)
           
            if(isempty(varargin))
                my_help('MRTQuant.RefitDKIWithMKCurve');
                return;
            end
            
            json.CallFunction = 'MRTQuant.RefitDKIWithMKCurve';
            json.Description = my_help('MRTQuant.RefitDKIWithMKCurve');
                        
            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end
            file_out = GiveValueForName(coptions,'output');
            if(isempty(file_out))
                error('Need to specify the output .mat file');
            end
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Quantification';
           
           mkcurve_fit_mat(file_in,file_out); 
           NiftiIO_basic.WriteJSONDescription('output',file_out(1:end-4),'props',json);
        end
        
        function PerformDTIBased_FiberTracking(varargin)
            % This function is deprecated. Please use MRTTrack.PerformDTIBased_FiberTracking
            warning('% This function is deprecated. Please use MRTTrack.PerformDTIBased_FiberTracking');            
        end
        
        function PerformCSD(varargin)
           % This function is deprecated. Please use MRTTrack.PerformCSD
            warning('% This function is deprecated. Please use MRTTrack.PerformCSD');            
        end   

        function PerformFODBased_FiberTracking(varargin)
           % This function is deprecated. Please use MRTTrack.PerformFODBased_FiberTracking
            warning('% This function is deprecated. Please use MRTTrack.PerformFODBased_FiberTracking');            
        end        

        function Perform_mFOD_FiberTracking(varargin)
           % This function is deprecated. Please use MRTTrack.Perform_mFOD_FiberTracking
            warning('% This function is deprecated. Please use MRTTrack.Perform_mFOD_FiberTracking');            
        end
        
        function MatMetrics2Nii(mat_file_in,dki_export)
        % This function exports the DTI/DKI metrics contained in a .mat
        % file to Nifti format. This will export fractional anisotropy (FA),
        % mean diffusivity (MD), axial diffusivity (L1), radial diffusivity
        % (RD). If the kurtosis export is enabled, mean kurtosis (MK),
        % axial kurtosis (AK), radial kurtosis (RK) and kurtosis anisotropy
        % (KA) will also be exported. Input arguments are:
        % mat_file_in: the target .mat file in ExploreDTI-like format
        % dki_export: whether to export kurtosis metrics, as 0 (default) or 1
            
            json.CallFunction = 'MRTQuant.MatMetrics2Nii';
            json.Description = my_help('MRTQuant.MatMetrics2Nii');

            list2export = [1 2 3 6 13 19];
            if(nargin > 1 && dki_export == 1)
                list2export = [list2export 46 47 48 49];
            end
            CLL = EDTI_Library.E_DTI_Complete_List_var;
            [fp,fn] = fileparts(mat_file_in);
            if(isempty(fp))
                fp = pwd;
            end
            KT = load(mat_file_in,'KT');
            has_dki = 0;
            if(isfield(KT,'KT'))
                has_dki = 1;
            end
            
            DTI_extensions = {'FA','MD','L1','RD','FE','FEFA'};
            DKI_extensions = {'MK','AK','RK','KA'};
            
            EDTI_Library.E_DTI_Convert_mat_2_nii(mat_file_in,fp,CLL(list2export));
            
            json.ReferenceFile = mat_file_in;
            json.ProcessingType = 'Metric';
            
            for metric_id=1:length(DTI_extensions)
                json.MetricDescription = CLL{metric_id};
                NiftiIO_basic.WriteJSONDescription('output',fullfile(fp,[fn '_' DTI_extensions{metric_id}]),'props',json);
            end
            if(has_dki == 1 && nargin > 1 && dki_export == 1)
                for metric_id=1:length(DKI_extensions)
                    json.MetricDescription = CLL{metric_id+length(DTI_extensions)};
                    NiftiIO_basic.WriteJSONDescription('output',fullfile(fp,[fn '_' DKI_extensions{metric_id}]),'props',json);
                end                
            end
        end
        
        function EnforceJSON(true_or_false)
        % Set a global flag which will force all functions to save
        % a descriptive .JSON file associated to the processing.
            global MRIToolkit;
            MRIToolkit.EnforceJSON = true_or_false;            
        end
        
        function EnforceNiiGz(true_or_false)
        % Set a global flag which will force all functions to save
        % compressed niftis.
            global MRIToolkit;
            MRIToolkit.EnforceNiiGz = true_or_false;
        end
        
        function ReplaceInvalidKurtosisPoints(true_or_false)
            % Sets a global flag which enables the replace of clearly bad
            % points (kurtosis < 0 or kurtosis >> 10) in the DKI fit.
            global MRIToolkit;
            MRIToolkit.DKI_cleanup = true_or_false;
        end
        
        function mrt_data = EDTI_Data_2_MRIToolkit(varargin)
        % This function loads the content of an ExploreDTI-like .mat file
        % and converts it to the data format used in MRIToolkit. mrt_data
        % is a struct with fields "img" (the 4D matrix), "bvals" (the b-values),
        % "bvecs" (the diffusion gradients), % "mask" (a mask derived from
        % EDTI). Input arguments:
        % mat_file: the .mat file to import

            if(isempty(varargin))
                my_help('MRTQuant.EDTI_Data_2_MRIToolkit');
                return;
            end  
            
            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end
            
            load(file_in,'DWI','NrB0','b','g','FA','VDims','eigval'); %,'bval','MDims','DT');
            
            B0s = EDTI_Library.E_DTI_Clean_up_B0s_2(DWI, ~isnan(FA), NrB0); % this part cleans the B0s where DWIs intensity> B0s intensity
            DWI(1:NrB0) = B0s;
            clear B0s;
            
            DWI = EDTI_Library.E_DTI_DWI_cell2mat(DWI); % all the volumes, B0s + DWIs
            mask = ~isnan(FA);
            bvals = round(sum(b(:,[1 4 6]),2)/100)*100;
            bvecs = [zeros(NrB0,3);g];
            
            mrt_data.img = single(DWI);
            mrt_data.bvals = bvals;
            mrt_data.bvecs = bvecs;
            mrt_data.mask = mask;
            mrt_data.VD = VDims;
        end        
        
        function EDTI_Data_2_Nii(varargin)
        % This function loads the content of an ExploreDTI-like .mat file
        % and converts it to .nii/.bval/.bvec 
        % Input arguments:
        % mat_file: the .mat file to import
        % output: the output prefix (adds .nii/.bval/.bvec)
            if(isempty(varargin))
                my_help('MRTQuant.EDTI_Data_2_MRIToolkit');
                return;
            end  
            
            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end
            
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                error('Need to specify the output prefix');
            end
            
            mrt_data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',file_in);
            MRTQuant.WriteNifti(mrt_data,[output '.nii']);
            
            bval = mrt_data.bvals';
            bvec = mrt_data.bvecs';
            
            fb = fopen([output '.bval'],'wt');
            for ij=1:length(bval)
               fprintf(fb,'%f ', bval(ij)); 
            end
            fclose(fb);
            
            fb = fopen([output '.bvec'],'wt');
            for ik=1:size(bvec,1)
                for ij=1:size(bvec,2)
                   fprintf(fb,'%f ',bvec(ik,ij)); 
                end
                fprintf(fb,'%s',newline);
            end
            fclose(fb);
        
        end

        function PerformMocoEPI(varargin)
        % This function performs the motion/eddy currents/EPI distortions
        % correction using Elastix. Input arguments:
        % mat_file: The ExpoloreDTI-like .mat file to correct
        % do_moco: 0-1 disabling/enabling the motion/eddy currents correction part
        % epi_tgt: The T1/T2 to use for EPI correction (via registration).
        %       If set, this automatically enables the registration step
        % constraint_epi: a vector in the form [0 1 0] enabling or disabling
        %       LR/AP/FH deformations
        % epi_reg_mode: Which image to use for registration toward the epi_tgt
        %       'b0' (use the first b=0s/mm2, default), 'avg_dwis' (use the average
        %       of all DWIs), 'fa' (use the fractional anisotropy)
        % epi_reg_type: 'linear' (rigid transformation) or 'nonlinear'
        % (default)
        % use_normcorr: 0-1 use normalized correlation in place of mutual
        %       information. The first is useful for images with similar
        %       contrasts (e.g. T1-FA), whereas mutual information typically works well
        %       both with similar and with different contrasts.
        % fit_mode: force a new fit_mode, ols, wls, nls, rekindle
        % output: output folder
        
            global MRIToolkit;         

            if(isempty(varargin))
                my_help('MRTQuant.PerformMocoEPI');
                return;
            end  

            json.CallFunction = 'MRTQuant.PerformMocoEPI';
            json.Description = my_help('MRTQuant.PerformMocoEPI');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';
            
            basic_info = load(file_in,'par');
            if(~isfield(basic_info,'par'))
%                 warning('This MAT file is missing essential information. Please, reprocess it with MRIToolkit');
                basic_info.par.clean_up_PIS = 1;
                basic_info.par.ROBUST_option = 1;
                basic_info.par.RE.rel_convergence = 1e-3;
                basic_info.par.RE.max_iter = 20;
                basic_info.par.RE.kappa = 6;
                basic_info.par.DKI_constraints.do_it = 1;
                basic_info.par.DKI_constraints.constr1 = [-Inf Inf];
                basic_info.par.DKI_constraints.constr2 = [-Inf Inf];
                basic_info.par.TE = 1;
            end
            
            option_value = GiveValueForName(coptions,'fit_mode');
            if(~isempty(option_value))
                fit_mode = option_value;
                if(strcmpi(fit_mode,'ols'))
                    basic_info.par.TE = 1;
                elseif(strcmpi(fit_mode,'wls'))
                    basic_info.par.TE = 2;
                elseif(strcmpi(fit_mode,'nls'))
                    basic_info.par.TE = 3;
                elseif(strcmpi(fit_mode,'rekindle'))
                    basic_info.par.TE = 4;
                end                    
            end
            
            par = basic_info.par;
            par.cust_mask.NS = '';
            par.cust_mask.TS = '';
            try
                mpar = load(file_in,'Mask_par');
                par.mask_P.NS.mfs = 5;
                par.mask_P.NS.NDWI = mpar.Mask_par.tune_NDWI;
                par.mask_P.NS.DWI = mpar.Mask_par.tune_DWI;
                par.mask_P.TS.mfs = 5;
                par.mask_P.TS.NDWI = mpar.Mask_par.tune_NDWI;
                par.mask_P.TS.DWI = mpar.Mask_par.tune_DWI;
            catch
                par.mask_P.NS.mfs = 5;
                par.mask_P.NS.NDWI = 0.7;
                par.mask_P.NS.DWI = 0.7;
                par.mask_P.TS.mfs = 5;
                par.mask_P.TS.NDWI = 0.7;
                par.mask_P.TS.DWI = 0.7;
            end
            par.R2D.type = 0;
            par.R2D.FN = '';
            par.R2D.contrast = 1;
            par.Num_iter = 1000;
            par.Num_samp = 2000;
            par.Hist_bin = 64;
            par.Num_Resol = 1;
            par.use_f_mask = 1;
            par.Regul = 0;
            par.suc = 1;
            par.suff.NS = '_MD_C_native.mat';
            par.suff.TS = '_MD_C_trafo.mat';
            par.DTI_f_in = file_in;
            par.no_GUI = 0;
            par.clean_up_PIS = 1;
            par.EPI.Num_iter = 1000;
            par.EPI.Num_samp = 20000;
            par.EPI.Hist_bin = 64;
            par.EPI.Num_Resol = 4;
            par.EPI.Grid_Spacing = [30 30 30];
            par.EPI.Deriv_Scales = [1 1 1];
            par.EPI.use_f_mask = 1;
            te = par.TE;
            par = rmfield(par,'TE');
            par.TE.NS = te;
            par.TE.TS = te;
            
            option_value = GiveValueForName(coptions,'do_moco');
            if(~isempty(option_value))
                if(option_value == 0)
                    par.DOF = 1;
                else
                    par.DOF = 3;
                end
            else
                par.DOF = 3;
            end
            
            option_value = GiveValueForName(coptions,'epi_tgt');
            if(~isempty(option_value))
                par.R2D.FN = option_value;
            else
                par.R2D.type = 0;
            end
            
            option_value = GiveValueForName(coptions,'epi_reg_type');
            if(~isempty(par.R2D.FN) && isempty(option_value))
                par.R2D.type = 3;
            else
                if(strcmpi(option_value,'linear'))
                    par.R2D.type = 1;
                elseif(strcmpi(option_value,'nonlinear'))
                    par.R2D.type = 3;
                end
            end
            
            option_value = GiveValueForName(coptions,'constraint_epi');
            if(~isempty(option_value))
                if(~ischar(option_value))
                    par.EPI.Deriv_Scales = [(option_value(1)) (option_value(2)) (option_value(3))];
                else
                    par.EPI.Deriv_Scales = [str2double(option_value(1)) str2double(option_value(2)) str2double(option_value(3))];
                end
            end
            
            option_value = GiveValueForName(coptions,'epi_reg_mode');
            if(~isempty(option_value))
                if(strcmpi(option_value,'b0'))
                    par.R2D.contrast = 1;
                elseif(strcmpi(option_value,'avg_dwis'))
                    par.R2D.contrast = 3;
                elseif(strcmpi(option_value,'fa'))
                    par.R2D.contrast = 2;
                else
                    error(['Unexpected epi_reg_mode value: ' option_value]);
                end
            end
            
            option_value = GiveValueForName(coptions,'output');
            if(~isempty(option_value))
                out_folder = option_value;
            else
                out_folder = fileparts(file_in);
                if(isempty(out_folder) || strcmp(out_folder,filesep))
                    out_folder = pwd;
                end
            end
            
            par.mask_P.TE.NS = par.TE;
            par.mask_P.TE.TS = par.TE;
            
            option_value = GiveValueForName(coptions,'use_normcorr');
            if(~isempty(option_value) && ...
                    ((isnumeric(option_value) && option_value == 1) ...
                    || isstring(option_value) && strcmpi(option_value,'1')))
                par.use_NC = 1;
            end
            
            par.temp_folder = fullfile(tempdir,['MRT_' num2str(randi(150)) '_' num2str(randi(150)) '_' num2str(randi(150))]);
            if(strcmp(par.temp_folder(end),'/'))
                par.temp_folder = par.temp_folder(1:end-1);
            end
            
            par.E_path = MRIToolkit.Elastix.ElastixCMD;
            par.T_path = MRIToolkit.Elastix.TransformixCMD;
            
            par.MDC_constr_fac = 6;
            
            par.Interpol = 1;
            
            if isempty(par.E_path) || isempty(par.T_path)
                if(~isdeployed)
                    par.suc=0;
                    disp('Cannot find an appropriate Elastix build...')
                    return;
                end
            end
            
            par.out_folder = out_folder;            
            
            % h_w = EDTI_Library.my_waitbar(0,'Checking data...');pause(0.5)
            DTI_files = EDTI_Library.E_DTI_SMECEPI_Get_input_DTI(par);
            % EDTI_Library.my_waitbar(1);close(h_w);pause(0.01);
            
            if isempty(DTI_files)
                disp('Could not find relevant *.mat files... ')
                return;
            end
            
            suc = EDTI_Library.E_DTI_do_initial_check_reg_tra_files(par.E_path,par.T_path);
            
            if suc==0
                if par.no_GUI==0
%                    EDTI_Library.my_msgbox('See command prompt for more info...','Error...','Modal')
                end
                return;
            end
            
            if(isfield(MRIToolkit,'EnforceNiiGz'))
                enforce_niigz = MRIToolkit.EnforceNiiGz;
            else
                enforce_niigz = false;
            end
            MRTQuant.EnforceNiiGz(false);
            
            TS = tic;
            for i=1:length(DTI_files)
                disp('Processing file:')
                disp([DTI_files{i} ' ...'])
                EDTI_Library.E_DTI_SMECEPI_Single(DTI_files{i},par)
                disp('Done!')
            end
            
            MRTQuant.EnforceNiiGz(enforce_niigz);

            json.parameters = par;
            [~,filename] = fileparts(file_in);
            out_name = fullfile(out_folder,[filename(1:end-4) '_MD_C_trafo']);
            NiftiIO_basic.WriteJSONDescription('output',out_name,'props',json);

            disp('Processing finished!')
            
            t=toc(TS);
            if t<3600
                m = t/60;
                disp(['Total computation time was ' num2str(m) ' minutes.'])
            elseif t>3600 && t<(3600*24)
                h = t/3600;
                disp(['Total computation time was ' num2str(h) ' hours.'])
            else
                d = t/(3600*24);
                disp(['Total computation time was ' num2str(d) ' days.'])
            end
        end
        
        function SelectVolumesWithBvals(varargin)
        % Select a subset of volumes in 4D .nii file based on their diffusion weighting. 
        % The function expects a .mat file, or a .bval/.bvec files couple or a .txt (b-matrix) file in the folder
        % of the .nii file and with the same name
        % Input arguments:
        % nii_file: the .nii file of the original data
        % mat_file: the .mat file of the original data
        % bvals: an array containing the b-values to select, e.g. [0 1000]
        % output: the new .nii file containing a subset of the volumes.
        % tol: the tolerance around the specified b-values
            if(isempty(varargin))
                my_help('MRTQuant.SelectVolumesWithBvals');
                return;
            end
        
            json.CallFunction = 'MRTQuant.SelectVolumesWithBvals';
            json.Description = my_help('MRTQuant.SelectVolumesWithBvals');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                file_in = GiveValueForName(coptions,'nii_file');
                if(isempty(file_in))
                    error('Need to specify the input .nii or .mat file');
                end       
            end
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';
            
            bvals_list = GiveValueForName(coptions,'bvals');
            if(isempty(bvals_list))
                error('Need to specify the bvals to keep');
            end
            json.bvals_list = bvals_list;
            
            outname = GiveValueForName(coptions,'output');
            if(isempty(outname))
                error('Need to specify the output .nii file');
            end
            tol = GiveValueForName(coptions,'tol');

            use_txt = 0;
            if(contains(file_in(end-4:end),'mat'))
               dta = load(file_in,'DWI','b','VDims');
               [bval,bvec] = MRTQuant.bval_bvec_from_b_Matrix(dta.b);
               bval = bval';
               bvec = bvec';
               data = EDTI_Library.E_DTI_DWI_cell2mat(dta.DWI);
               VD = dta.VDims;
            else
                [data,VD] = EDTI_Library.E_DTI_read_nifti_file(file_in);
                if(exist([file_in(1:end-4) '.bval'],'file'))
                    bval = load([file_in(1:end-4) '.bval']);
                    bvec = load([file_in(1:end-4) '.bvec']);
                else
                    use_txt = 1;
                    bmat = load([file_in(1:end-4) '.txt']);
                    [bval,~] = MRTQuant.bval_bvec_from_b_Matrix(bmat);
                end
            end
            IX = false(size(bval));
            for ij=1:length(bvals_list)
                if(~isempty(tol))
                    u_tol = tol;
                else
                    u_tol = max(1,0.1*bvals_list(ij));
                end
                IX(abs(round(bval)-bvals_list(ij)) <= u_tol) = true;
            end     
            data = data(:,:,:,IX);
            bval = bval(IX);
            if(use_txt == 1)
                bmat = bmat(IX,:);
                save([outname(1:end-4) '.txt'],'bmat','-ascii');
            else
                bvec = bvec(:,IX);
                save([outname(1:end-4) '.bvec'],'bvec','-ascii');
                save([outname(1:end-4) '.bval'],'bval','-ascii');
            end
            EDTI_Library.E_DTI_write_nifti_file(data,VD,outname);

            NiftiIO_basic.WriteJSONDescription('output',outname(1:end-4),'props',json);
        end
        
        function SelectVolumesWithIndices(varargin)
        % Select a subset of volumes in 4D .nii file based on their 4d
        % order.
        % The function expects a .bval/.bvec files couple or a .txt (b-matrix) file in the folder
        % of the .nii file and with the same name 
        % Input arguments:
        % nii_file: the .nii file of the original data
        % indices: an array containing the b-values to select, e.g. [1 3 7 15 ...]
        % output: the new .nii file containing a subset of the volumes.
            if(isempty(varargin))
                my_help('MRTQuant.SelectVolumesWithIndices');
                return;
            end
        
            json.CallFunction = 'MRTQuant.SelectVolumesWithIndices';
            json.Description = my_help('MRTQuant.SelectVolumesWithIndices');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Need to specify the input .nii file');
            end       
            
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';
            
           indices_list = GiveValueForName(coptions,'indices');
            if(isempty(indices_list))
                error('Need to specify the indices to select');
            end
            json.indices_list = indices_list;
            
            outname = GiveValueForName(coptions,'output');
            if(isempty(outname))
                error('Need to specify the output .nii file');
            end

            use_txt = 0;
            if(exist([file_in(1:end-4) '.bval'],'file'))
                bval = load([file_in(1:end-4) '.bval']);
                bvec = load([file_in(1:end-4) '.bvec']);
            else
                use_txt = 1;
                bmat = load([file_in(1:end-4) '.txt']);
                [bval,~] = MRTQuant.bval_bvec_from_b_Matrix(bmat);
            end
            IX = indices_list;   
            [data,VD] = EDTI_Library.E_DTI_read_nifti_file(file_in);
            data = data(:,:,:,IX);
            bval = bval(IX);
            if(use_txt == 1)
                bmat = bmat(IX,:);
                save([outname(1:end-4) '.txt'],'bmat','-ascii');
            else
                bvec = bvec(:,IX);
                save([outname(1:end-4) '.bvec'],'bvec','-ascii');
                save([outname(1:end-4) '.bval'],'bval','-ascii');
            end
            EDTI_Library.E_DTI_write_nifti_file(data,VD,outname);

            NiftiIO_basic.WriteJSONDescription('output',outname(1:end-4),'props',json);
        end
        
        function [nrb0s,nrdwis,b0s_indices,dwis_indices] = GetNumOfB0sDWIs(varargin)
        % Gets the number of b=0s/mm2 and of dwis from a .mat or .nii
        % Input arguments (nii_file or mat_file are exclusive):
        % nii_file: .nii to inspect. a .txt or .bval / .bvec in the same folder and
        %       with the same name specifying the diffusion gradients are expected.
        % mat_file: .mat file to inspect.
        % b0_tol: the minimum b=0s/mm2 to be treated as 0. default = 5s/mm2
            if(isempty(varargin))
                my_help('MRTQuant.GetNumOfB0sDWIs');
                nrb0s = [];
                nrdwis = [];
                b0s_indices = [];
                dwis_indices = [];
                return;
            end

            coptions = varargin;
            b0_tol = GiveValueForName(coptions,'b0_tol');
            if(isempty(b0_tol))
                b0_tol = 5;
            end
            
            file_in = GiveValueForName(coptions,'nii_file');
            if(~isempty(file_in))
                if(exist([file_in(1:end-4) '.bval'],'file'))
                    bvals = load([file_in(1:end-4) '.bval']);
                else
                    bvals = load([file_in(1:end-4) '.txt']);
                    bvals = sum(bvals(:,[1 4 6]),2);
                end
                b0s_indices = find(bvals <= b0_tol);
                dwis_indices = find(bvals > b0_tol);
                nrb0s = length(b0s_indices);
                nrdwis = length(dwis_indices);
                return;
            end        
            
            file_in = GiveValueForName(coptions,'mat_file');
            if(~isempty(file_in))
                bvals = load(file_in(1:end-4), 'b');
                bvals = sum(bvals.b(:,[1 4 6]),2);
                b0s_indices = find(bvals <= b0_tol);
                dwis_indices = find(bvals > b0_tol);
                nrb0s = length(b0s_indices);
                nrdwis = length(dwis_indices);
                return;
            end   
            
        end
       
        function PerformGibbsRingingCorrection(varargin)
        % Performs the Gibbs ringing correction on b=0s/mm2 images with a
        % TV filter. Input arguments:
        % nii_file: .nii file to correct (sorted in ascending b-value
        % order)
        % output: Name of the output .nii
            if(isempty(varargin))
                my_help('MRTQuant.PerformGibbsRingingCorrection');
                return;
            end

            json.CallFunction = 'MRTQuant.PerformGibbsRingingCorrection';
            json.Description = my_help('MRTQuant.PerformGibbsRingingCorrection');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Missing mandatory argument nii_file');
            end
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';
            
            output = GiveValueForName(coptions,'output');
            if(isempty(file_in))
                error('Missing mandatory argument output');
            end
            EDTI_Library.GibbsRingingCorrection(file_in,output);
            
            if(exist([file_in(1:end-4) '.txt'],'file'))
                copyfile([file_in(1:end-4) '.txt'],[output(1:end-4) '.txt']);
            end
            
            NiftiIO_basic.WriteJSONDescription('output',output(1:end-4),'props',json);
        end
        
        function SortNiiWRTbval(varargin)
            % Sort a .nii according to the diffusion weighting in ascending
            % order. Input arguments:
            % nii_file: the .nii file to sort
            % output: the output .nii file
            % txt_file: the b-matrix file (optional). If not specified, the
            %   function will look for a .txt or .bval/.bvec with the same name
            %   of the .nii     
            if(isempty(varargin))
                my_help('MRTQuant.SortNiiWRTbval');
                return;
            end            
            
            json.CallFunction = 'MRTQuant.SortNiiWRTbval';
            json.Description = my_help('MRTQuant.SortNiiWRTbval');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Missing mandatory argument nii_file');
            end
            json.ReferenceFile = file_in;
            json.ProcessingType = 'Preprocessing';
            
            bmat_file = GiveValueForName(coptions,'txt_file');
            if(isempty(bmat_file))
                bmat_file = [file_in(1:end-4) '.txt'];
            end
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                output = [file_in(1:end-4) '_sorted.nii'];
            end            
            
            data = MRTQuant.LoadNifti(file_in);
            use_bmat = 0;
            if(exist(bmat_file,'file'))
                use_bmat = 1;
            end
            
            if(use_bmat == 1)
                bmat = load(bmat_file);
                bvals = sum(bmat(:,[1 4 6]),2);
                [~,IX] = sort(bvals,'ascend');
                bmat = bmat(IX,:);
                save([output(1:end-4) '.txt'],'bmat','-ascii');
            else
                bvals = load([file_in(1:end-4) '.bval']);
                bvecs = load([file_in(1:end-4) '.bvec']);
                [bvals,IX] = sort(bvals,'ascend');
                bvecs = bvecs(:,IX);
                save([output(1:end-4) '.bval'],'bvals','-ascii');
                save([output(1:end-4) '.bvec'],'bvecs','-ascii');                
            end
            
            data.img = data.img(:,:,:,IX);
            MRTQuant.WriteNifti(data,output);
            
            NiftiIO_basic.WriteJSONDescription('output',output(1:end-4),'props',json);
        end
        
        function DiscardEmptyVolumes(varargin)
        % Checks a 4D .nii for empty volumes and discards them
        % arguments:
        % nii_file: the (d)MRI .nii file
        % bmat: the corresponding b-matrix in .txt
        % output: the output .nii file
        
            if(isempty(varargin))
                my_help('MRTQuant.DiscardEmptyVolumes');
                return;
            end 
            
            json.CallFunction = 'MRTQuant.DiscardEmptyVolumes';
            json.Description = my_help('MRTQuant.DiscardEmptyVolumes');
            
            coptions = varargin;
            nii_file = GiveValueForName(coptions,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            json.ReferenceFile = nii_file;
            json.ProcessingType = 'Preprocessing';
            
            bmat = GiveValueForName(coptions,'bmat');
            if(~isempty(bmat))
                bmat = load(bmat);
            end
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end

            data = MRTQuant.LoadNifti(nii_file);
            [sx,sy,sz,st] = size(data.img);
            data.img = reshape(data.img,sx*sy*sz,st);
            V = sum(data.img);
            BI = V == 0;
            data.img = reshape(data.img,sx,sy,sz,st);
            if(sum(BI) ~= 0)
                data.img(:,:,:,BI == 1) = [];
                if(~isempty(bmat))
                   bmat(BI == 1,:) = [];
                   save([output(1:end-4) '.txt'],'bmat','-ascii');
                end
                MRTQuant.WriteNifti(data,output);
            end
            
            NiftiIO_basic.WriteJSONDescription('output',output(1:end-4),'props',json);
        end
        
        function res = AverageNormalizedResiduals(varargin)
        % Compute voxel-wise residuals using the DTI/DKI model. Input
        % arguments:
        % mat_file: the ExploreDTI-like .mat file to compute residuals from
        % normalize: 0-1 divide all signals by the correspondent b = 0s/mm2
        % output: if set, save the residuals to the specified .nii file
            
            % ADD JSON HANDLING
        
            if(isempty(varargin))
                my_help('MRTQuant.AverageNormalizedResiduals');
                return;
            end 
        
            coptions = varargin;
            mat_file = GiveValueForName(coptions,'mat_file');
            if(isempty(mat_file))
                error('Missing mandatory argument mat_file');
            end
            normalize = GiveValueForName(coptions,'normalize');
            if(isempty(normalize))
                normalize = 1;
            end
            load(mat_file,'DWI','DWIB0','DT','b','VDims','KT');
            DWI = single(E_DTI_DWI_cell2mat(DWI));
            if(normalize == 1)
                for iz=1:size(DWI,4)
                    sl = DWI(:,:,:,iz)./(DWIB0+eps);
                    sl(DWIB0 == 0) = 0;
                    DWI(:,:,:,iz) = sl;
                end
            end
            DT = EDTI_Library.E_DTI_DWI_cell2mat(DT);
            siz = size(DT);
            if(exist('KT','var') < 1 || isempty(KT))
                DT = reshape(DT,siz(1)*siz(2)*siz(3),size(DT,4));
                pred = exp(DT*(-b'));
            end
            pred = reshape(pred,[siz(1:3) size(b,1)]);
            if(normalize == 0)
                for iz=1:size(pred,4)
                    pred(:,:,:,iz) = DWIB0.*pred(:,:,:,iz);
                end
            end
            res = mean(abs(DWI-pred),4);
            
            output = GiveValueForName(coptions,'output');
            if(~isempty(output))
                EDTI_Library.E_DTI_write_nifti_file(res,VDims,output);
            end
            
        end
                                
        function temp_mat_location = QuickNiiBvalBvecToMat(varargin)
            % Create a .MAT from .nii and .bval/.bvec for internal processing
            % in a temporary directory. Input arguments:
            % nii_file: the target .nii file
            % bval_file: the compation .bval fie
            % bvec_file: the companion .bvec file
            % grad_perm: optional, how to permute the diffusion gradients
            % grad_flip: optional, how to flip the sign of the gradients
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
            
            MRTQuant.b_Matrix_from_bval_bvec('bval_file',bval_file,'bvec_file',bvec_file,...
                'output',[dest_basename '.txt']);
            
            if(~isempty(perm) && ~isempty(flip))
                MRTQuant.PerformDTI_DKIFit('nii_file',file_in,'txt_file',[dest_basename '.txt'],...
                    'output',[dest_basename '.mat'],'grad_perm',perm,'grad_flip',flip);
            else
                MRTQuant.PerformDTI_DKIFit('nii_file',file_in,'txt_file',[dest_basename '.txt'],...
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
            vol = MRTQuant.LoadNifti(data_name);
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
                MRTQuant.WriteNifti(OUT,[save_prefix '_denoised.nii']);
                OUT.img = noise_map;
                MRTQuant.WriteNifti(OUT,[save_prefix '_noisemap.nii']);
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
                my_help('MRTD.PerformSpectralDeconvolution');
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
            % fit_dt: (default 0), whether to fit a tensor for the tissue
            % b_thresh: The b-value after which IVIM is "negligible"
            % (typically b=200s/mm2)
            % output: the output prefix name (suffixes and .nii will be
            % added automatically)   
            % do_geoavg: 0-1, whether to perform geometric averaging
            % (default 1)
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = MRTQuant.LoadNifti(data_name);

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
                    [data.bvals,data.bvecs] = MRTQuant.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
            end
            
            mask_name = GiveValueForName(coptions,'mask_file');
            if(isempty(mask_name))
                error('Need to specify the mask file (mask_file)');
            end

            fit_dt = GiveValueForName(coptions,'fit_dt');
            if(isempty(fit_dt))
                fit_dt = 0;
            end
            
            do_geoavg = GiveValueForName(coptions,'do_geoavg');
            if(isempty(do_geoavg))
                do_geoavg = 1;
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

            mask = MRTQuant.LoadNifti(mask_name);
            data.mask = mask.img;
            
            if(do_geoavg == 1)
                data = QuickGeoAverage(data);
            end
            
            if(fit_dt == 0)
                [Dhigh,Dlow,f] = SegmentedIVIMFit2Comp(data,parameters);
                MRTQuant.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
            else
                [DT,Dlow,f] = SegmentedIVIMFit2Comp_DT(data,parameters);
                [~, FA, ~, ~, eigval] = EDTI_Library.E_DTI_eigensystem_analytic(EDTI_Library.E_DTI_DWI_mat2cell(DT.img));
                out_nii = f;
                out_nii.img = FA;
                MRTQuant.WriteNifti(out_nii,[output_prefix '_FA.nii']);
                out_nii.img = mean(eigval,4);
                MRTQuant.WriteNifti(out_nii,[output_prefix '_MD.nii']);                
            end
            MRTQuant.WriteNifti(f,[output_prefix '_fivim.nii']);
            MRTQuant.WriteNifti(Dlow,[output_prefix '_Dstar.nii']);
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
            % correct_data: 0 or 1. Whether to remove the FW signal from the data
            % (which can be later fit with DTI to get FW-corrected FA/MD)
            % do_geoavg: 0-1, whether to perform geometric averaging
            % (default 1)
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = MRTQuant.LoadNifti(data_name);

                mask_name = GiveValueForName(coptions,'mask_file');
                if(isempty(mask_name))
                    error('Need to specify the mask file (mask_file)');
                end
                mask = MRTQuant.LoadNifti(mask_name);
                data.mask = mask.img;
            
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
                    [data.bvals,data.bvecs] = MRTQuant.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
                FA = load(mat_name,'FA');
                data.mask = ~isnan(FA.FA);
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

            correct_data = GiveValueForName(coptions,'correct_data');
            if(isempty(correct_data))
                correct_data = 0;
            end            
            
            do_geoavg = GiveValueForName(coptions,'do_geoavg');
            if(isempty(do_geoavg))
                do_geoavg = 1;
            end

            if(do_geoavg == 1)
                data = QuickGeoAverage(data);
            end
            
            [Dhigh,f,S0] = FWFit2Comp(data,parameters);

            MRTQuant.WriteNifti(f,[output_prefix '_ffw.nii']);
            MRTQuant.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
            
            if(correct_data == 1)
               for ij=1:size(data.img,4)
                   data.img(:,:,:,ij) = max(data.img(:,:,:,ij) - S0.img.*f.img.*exp(-3e-3*data.bvals(ij)),0); 
               end
               MRTQuant.WriteNifti(data,[output_prefix '_FWcorrected.nii']);
               bmat = MRTQuant.b_Matrix_from_bval_bvec('bval',data.bvals,'bvec',data.bvecs);
               save([output_prefix '_FWcorrected.txt'],'bmat','-ascii');
            end
            
        end                

        function PerformIVIM_Fit(varargin)
            % Perform an IVIM + ADC non-linear fit.
            % The input MUST have multiple diffusion weightings (2+) including low b-values (b<=200s/mm2)
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
            % do_geoavg: 0-1, whether to perform geometric averaging
            % (default 1)
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = MRTQuant.LoadNifti(data_name);

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
                    [data.bvals,data.bvecs] = MRTQuant.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
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
            
            do_geoavg = GiveValueForName(coptions,'do_geoavg');
            if(isempty(do_geoavg))
                do_geoavg = 1;
            end
            
            mask = MRTQuant.LoadNifti(mask_name);
            data.mask = mask.img;
            
            if(do_geoavg == 1)
                data = QuickGeoAverage(data);
            end
            
            [Dhigh,fivim,DStar] = IVIMFit2Comp(data,parameters);

            MRTQuant.WriteNifti(fivim,[output_prefix '_fivim.nii']);
            MRTQuant.WriteNifti(DStar,[output_prefix '_Dstar.nii']);
            MRTQuant.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
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
            % min_bval: Minimum b-value to use
            % max_bval: Maximum b-value to use
            % output: the output prefix name (suffixes and .nii will be
            % added automatically)   
            % do_geoavg: 0-1, whether to perform geometric averaging
            % (default 1)
            
            coptions = varargin;

            mat_name = GiveValueForName(coptions,'mat_file');           
            
            if(isempty(mat_name))
                data_name = GiveValueForName(coptions,'nii_file');
                if(isempty(data_name))
                    error('Need to specify the target file (nii_file)');
                end
                data = MRTQuant.LoadNifti(data_name);

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
                    [data.bvals,data.bvecs] = MRTQuant.bval_bvec_from_b_Matrix(gmat);
                end
            else
                data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',mat_name);
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
            
            do_geoavg = GiveValueForName(coptions,'do_geoavg');
            if(isempty(do_geoavg))
                do_geoavg = 1;
            end
            
            mask = MRTQuant.LoadNifti(mask_name);
            data.mask = mask.img;
            
            if(do_geoavg == 1)
                data = QuickGeoAverage(data);
            end
            
            [Dhigh,Dlow,ffw,fivim] = IVIMFWFit3Comp(data,parameters);

            MRTQuant.WriteNifti(ffw,[output_prefix '_ffw.nii']);
            MRTQuant.WriteNifti(fivim,[output_prefix '_fivim.nii']);
            MRTQuant.WriteNifti(Dlow,[output_prefix '_Dstar.nii']);
            MRTQuant.WriteNifti(Dhigh,[output_prefix '_Dhigh.nii']);
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

% Helper: validate all the input parameters
function value = ValidateInputParameters(coptions,supported_params)
for ij=1:2:length(coptions)
    found = -1;
    for ik=1:length(supported_params)
       if(strcmpi(coptions{ij},supported_params{ik}))
           found = ik;
           break
       end
    end
    if(found == -1)
       error(['Unrecognized input parameter ' coptions{ij}]); 
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
   data.img(data.img < 0 | isnan(data.img)) = 0;
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
    data.img = single(data.img);
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

% Segmented IVIM + tissue (isotropic) fit
function [DT,Dlow,f] = SegmentedIVIMFit2Comp_DT(data,parameters)
    IX_high = data.bvals >= parameters.bthresh & data.bvals <= parameters.max_bval;
    IX_low = data.bvals >= parameters.min_bval & data.bvals < parameters.bthresh;   
    
    siz = size(data.img);
    data.img = log(reshape(data.img,siz(1)*siz(2)*siz(3),siz(4)));
%     data.img = (reshape(data.img,siz(1)*siz(2)*siz(3),siz(4)));
    
    bmat = MRTQuant.b_Matrix_from_bval_bvec('bval',data.bvals,'bvec',data.bvecs);
    Xh = [-bmat(IX_high,:) ones(sum(IX_high==1),1)];
    DTp = Xh\data.img(:,IX_high)';
%     DTp = EDTI_Library.E_DTI_WLLS_WW(data.img(:,IX_high),Xh');
    DTp = permute(DTp,[2 1]);
    DTp(data.mask(:) == 0,:) = 0;
    
    DT.img = reshape(DTp(:,1:6),[siz(1:3) 6]);
    DT.VD = data.VD;

    Xl =[-data.bvals(IX_low) ones(sum(IX_low==1),1)];
    LowD = Xl\data.img(:,IX_low)';
    Dlow.img = reshape(LowD(1,:),siz(1:3));
    Dlow.VD = data.VD;
    Dlow.img(data.mask == 0) = 0;

    data.img = exp(data.img);
    f.img = abs(data.img(:,1) - exp(DTp(:,7)))./data.img(:,1);
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
function [Dhigh,f,S0] = FWFit2Comp(data,parameters)
    data.img = double(data.img);
    IX_high = data.bvals >= parameters.min_bval & data.bvals <= parameters.max_bval;

    x0 = [1 0.05 0.7e-3];
    lb = [0 0 0];
    ub = [1 1 2e-3];
    
    Dhigh.img = zeros(size(data.img(:,:,:,1)));
    Dhigh.VD = data.VD;
    f.img = zeros(size(data.img(:,:,:,1)));
    f.VD = data.VD;
    S0.img = zeros(size(data.img(:,:,:,1)));
    S0.VD = data.VD;
    
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
        S0.img(points2fit(pi)) = p(1);
    end
    toc
end

% IVIM + D model
function [out,S] = IVIM2CompartmentsModel(x0,params)
    S = x0(1)*((1-x0(2))*exp(-params.b*x0(3))+x0(2)*exp(-params.b*x0(4)));
    out = params.S - S;
end

% FW + Tissue (isotropic) non-linear fit
function [Dhigh,f,DStar,S0] = IVIMFit2Comp(data,parameters)
    data.img = double(data.img);
    IX_high = data.bvals >= parameters.min_bval & data.bvals <= parameters.max_bval;

    x0 = [1 0.05 0.7e-3 50e-3];
    lb = [0 0 0 5e-3];
    ub = [1 1 3e-3 300e-3];
    
    Dhigh.img = zeros(size(data.img(:,:,:,1)));
    Dhigh.VD = data.VD;
    f.img = zeros(size(data.img(:,:,:,1)));
    f.VD = data.VD;
    S0.img = zeros(size(data.img(:,:,:,1)));
    S0.VD = data.VD;
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
            disp(['IVIM fit: ' num2str(pi/length(points2fit)*100) '%']);
        end
        params.S = data.img(IX_high,points2fit(pi));
        x0(1) = params.S(1);
        ub(1) = 2*params.S(1);
        p = lsqnonlin(@IVIM2CompartmentsModel,x0,lb,ub,options,params);
        DStar.img(points2fit(pi)) = p(4);
        Dhigh.img(points2fit(pi)) = p(3);
        f.img(points2fit(pi)) = p(2);
        S0.img(points2fit(pi)) = p(1);
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
    data.img = double(data.img);
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

% This function implements the MK curve fit and is not part of the 
% original code of ExploreDTI. (10.1016/j.neuroimage.2019.04.015)
function mkcurve_fit_mat(input_file,output_file)
    data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',input_file);
    gradient_info = load(input_file,'b','g','NrB0');
    tensors = load(input_file,'KT','DT','FA');
%     if(~isfield(tensors,'outlier'))
%         tensors.outlier = false(size(data.img));
%     end
    siz = size(data.img);
    V = reshape(data.img,siz(1)*siz(2)*siz(3),siz(4));
    V(isnan(tensors.FA(:)),:) = 0;

    tic
    g = [zeros(gradient_info.NrB0,3);gradient_info.g];
    b = gradient_info.b;
    bval = sum(b(:,[1 4 6]),2);
    
    b_kurt = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];

    b_kurt = repmat((1/54)*(b(:,1)+b(:,4)+b(:,6)).^2,[1 15]).*b_kurt;

    b_final = [ones(size(data.img,4),1) -b b_kurt];
    Xd = pinv(b_final);

    K1 = Xd*log(V');
    DT = reshape(K1(2:7,:)',[siz(1:3) 6]);
    KT = reshape(K1(8:end,:)',[siz(1:3) 15]);
    MD = (DT(:,:,:,1)+DT(:,:,:,4)+DT(:,:,:,6))/3;
    KT = KT./((MD*3).^2);

    DT = EDTI_Library.E_DTI_DWI_mat2cell(DT);
    KT = EDTI_Library.E_DTI_DWI_mat2cell(KT);

%     MK = EDTI_Library.E_DTI_Mean_Kurtosis(KT,DT);

    A = [g(:,1).^4 g(:,2).^4 g(:,3).^4 ...
        4*(g(:,1).^3).*g(:,2) 4*(g(:,1).^3).*g(:,3) ...
        4*(g(:,2).^3).*g(:,1) 4*(g(:,2).^3).*g(:,3) ...
        4*(g(:,3).^3).*g(:,1) 4*(g(:,3).^3).*g(:,2) ...
        6*(g(:,1).^2).*(g(:,2).^2) ...
        6*(g(:,1).^2).*(g(:,3).^2) ...
        6*(g(:,2).^2).*(g(:,3).^2) ...
        12*g(:,2).*g(:,3).*(g(:,1).^2) ...
        12*g(:,1).*g(:,3).*(g(:,2).^2) ...
        12*g(:,1).*g(:,2).*(g(:,3).^2)];

    B = [g(:,1).^2 g(:,2).^2 g(:,3).^2 2.*g(:,1).*g(:,2) 2.*g(:,1).*g(:,3) 2.*g(:,2).*g(:,3)];

    [sx,sy,sz,st] = size(data.img);

%     MK_Q = zeros(size(data.img(:,:,:,1)));
%     FA_Q = zeros(size(data.img(:,:,:,1)));
%     MD_Q = zeros(size(data.img(:,:,:,1)));
%     KA_Q = zeros(size(data.img(:,:,:,1)));
    NEW_S0 = zeros(size(data.img(:,:,:,1)));
    DT_Q = EDTI_Library.E_DTI_DWI_cell2mat(DT);
    KT_Q = EDTI_Library.E_DTI_DWI_cell2mat(KT);

    DT_Q = reshape(DT_Q,sx*sy*sz,size(DT_Q,4));
    KT_Q = reshape(KT_Q,sx*sy*sz,size(KT_Q,4));
    
    lambda = 0.5;

    S0 = data.img(:,:,:,1);
    min_val = max(S0(:))*0.01;

    data.img = reshape(data.img,sx*sy*sz,st);
%     outliers = reshape(tensors.outlier,sx*sy*sz,st);
    
    parfor x=1:size(data.img,1)
        S = [];
        GP = [];
        So = squeeze(data.img(x,:))';
        if((So(1) < min_val) || any(~isfinite(So)) || isnan(tensors.FA(x)))
            continue
        end
        l_bval = bval;%(outliers(x,:)==0);
        s0_val = mean(So(l_bval <= 1));
        s0_range = linspace(0.5*s0_val,s0_val*2,100);
%         [~,C_IX] = min(abs(s0_range-s0_val));

        MK = zeros(length(s0_range),1);

        for s0_id=1:length(s0_range)

            S = So;
            S(l_bval <= 1) = s0_range(s0_id);

            GP = S' ~= 0;% & outliers(x,:) == 0;
            X = Xd(:,GP)*log(S(GP));

            DT = X(2:7,:);
            MD = (DT(1)+DT(4)+DT(6))/3;
            KT = X(8:end,:)./((MD*3).^2);

            dt = DT([1 4 6 2 3 5]);
            MDsq = MD^2;

            for i=gradient_info.NrB0+1:size(g,1)
                if(GP(i) == 0)
                    continue
                end
                
                dum_adc = (B(i,:)*dt).^2;
                dum_kt = A(i,:)*KT;

                delta = (MDsq.*dum_kt)./dum_adc;
                if(isnan(delta))
                    disp('Debug');
                    continue
                end

                MK(s0_id) = MK(s0_id) + delta;

            end

            MK(s0_id) = MK(s0_id)/(size(g,1)-sum(GP(l_bval > 1) == 0));
        end

        last_negative = find(MK<0,1,'last');
        MKc = MK;
        MKc(1:last_negative) = NaN;
        [~,maxmks] = max(MKc);

        zeromks = find(diff(sign(MK)) ~= 0);
        if(length(zeromks) > 1)
            [~,index] = min(abs(s0_range(zeromks)-s0_range(maxmks)));
            zeromks = zeromks(index);
        end

%         MK(1:zeromks) = NaN;

%         corr_mk = C_IX;
        corr_s0 = s0_val;

        if((~isempty(zeromks) && ~isempty(maxmks)) && zeromks < maxmks)
            corr_s0 = (1-lambda)*s0_range(zeromks)+lambda*s0_range(maxmks);
%             if(s0_val < corr_s0)
%                 [~,corr_mk] = min(abs(s0_range-corr_s0));
%             end
        end

        S(1:gradient_info.NrB0) = corr_s0;
        X = Xd(:,GP)*log(S(GP));

        DT = X(2:7,:);
        MD = (DT(1)+DT(4)+DT(6))/3;
        KT = X(8:end,:)./((MD*3).^2);

        KT_Q(x,:) = KT;
        DT_Q(x,:) = DT;

%         dt = DT([1 4 6 2 3 5]);
%         MDsq = MD^2;

%         for i=gradient_info.NrB0+1:size(g,1)
% 
%             dum_adc = (B(i,:)*dt).^2;
%             dum_kt = A(i,:)*KT;
%             dum = (MDsq.*dum_kt)./dum_adc(:);

%             MK_Q(x) = MK_Q(x) + dum;
%             KA_Q(x) = KA_Q(x) + dum^2;
%         end
%         MK_Q(x) = MK_Q(x)/size(g,1);    
%         KA_Q(x) = KA_Q(x)/size(g,1);    

%         DT = [DT(1) DT(2)/2 DT(3)/2;
%               DT(2)/2 DT(4) DT(5)/2;
%               DT(3)/2 DT(5)/2 DT(6)];
%         [EVEC,EVAL] = eig(DT);
%         EVAL = diag(EVAL);
%         FA_Q(x) = sqrt(1.5*sum((EVAL-mean(EVAL)).^2)./sum(EVAL.^2));
%         MD_Q(x) = MD;

%         if(isnan(MK(corr_mk)))
%             MK_Q(x) = MK(end);
%         end
%         KA_Q(x) = sqrt(abs(KA_Q(x)-MK_Q(x)^2));
        NEW_S0(x) = corr_s0;
    end
    toc

    DT_Q = reshape(DT_Q,sx,sy,sz,size(DT_Q,2));
    KT_Q = reshape(KT_Q,sx,sy,sz,size(KT_Q,2));
    
    DT = EDTI_Library.E_DTI_DWI_mat2cell(DT_Q);
    KT = EDTI_Library.E_DTI_DWI_mat2cell(KT_Q);
    DWIB0 = NEW_S0;
    try
        copyfile(input_file,output_file);
    catch
    end
    save(output_file,'DT','KT','DWIB0','-append');

end

function out = my_help(fname)
    if(isdeployed)
        out = fname;
    else
        out = help(fname);
    end
end
