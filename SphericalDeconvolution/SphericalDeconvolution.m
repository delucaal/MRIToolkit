%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



% Spherical deconvolution class - including GRL and mFOD
% A. De Luca - UMC Utrecht - 28/10/2019 - alberto@isi.uu.nl -
% alberto.deluca.06@gmail.com
% First version: 28/10/2019

classdef SphericalDeconvolution < handle
    properties
        data;
        n_isotropic;
        n_anisotropic;
        deconv_method;
        inner_shells_weight;
        rf_indices_lr;
        rf_indices_hr;
        rf_models;
        L2LSQ_reg;
        LRKernel;
        HRKernel;
        nDirections;
        NN_L, NN_H;
        data_size;
    end
    
    methods
        
        function obj = SphericalDeconvolution(varargin)
            % Class constructor. Accepts only 1 optional argument:
            % data, a structure with fields "img" (the 4D data matrix), "bvals"
            % the b-values, "bvecs" the gradient vectors
            
            obj.data = [];
            obj.n_anisotropic = 0;
            obj.n_isotropic = 0;
            obj.deconv_method = [];
            obj.inner_shells_weight = 0.2;
            obj.rf_indices_lr = {};
            obj.rf_indices_hr = {};
            obj.LRKernel = [];
            obj.HRKernel = [];
            obj.L2LSQ_reg = 0.1;
            obj.nDirections = 300;
            obj.NN_L = 0;
            obj.NN_H = 0;
            obj.rf_models = {};
            
            InputPairs = ParseInputKeys(varargin);
            for ij=1:size(InputPairs,1)
                if(strcmpi(InputPairs{ij,1},'data'))
                    obj.data = InputPairs{ij,2};
                    obj.data_size = size(obj.data.img);
                end
            end
        end
        
        function setInnerShellWeighting(obj,weight)
            %       % How the inner shells are weighted in the deconvolution (0-1)
            obj.inner_shells_weight = weight;
        end
        
        function boolean = isInitialized(obj)
            % Returns true if everything is set properly and ready to go
            boolean = (obj.n_anisotropic + obj.n_isotropic ~= 0) & ...
                SphericalDeconvolution.isSupportedMethod(obj.deconv_method) & ~isempty(obj.data) & ...
                ~isempty(obj.LRKernel) & ~isempty(obj.HRKernel) & obj.nDirections > 0;
            if(strcmpi(obj.deconv_method,'dRL'))
                boolean = boolean & obj.NN_L ~= 0 & obj.NN_H ~= 0;
            end
        end
        
        function obj = AddIsotropicRF(obj,D)
            % Add an isotropic response function (i.e. CSF) to the
            % deconvolution matrix. These should be add only after all
            % anisotropic RFs have been added.
            
            if(isempty(obj.data))
                warning('Cannot initialize RFs if data has not been set. Bailing out.');
                return
            end
            obj.n_isotropic = obj.n_isotropic+1;
            IndexHR = size(obj.HRKernel,2)+1;
            IndexLR = size(obj.LRKernel,2)+1;
            [~,lLRKernel,lHRKernel] = mDRLMT_MakeDKIKernel_multicomp(obj.data,10,[1.7e-3 0.5e-3 0.5e-3],0,D,1);
            if(~isempty(obj.LRKernel))
                obj.LRKernel(:,IndexLR) = lLRKernel(:,end);
                obj.HRKernel(:,IndexHR) = lHRKernel(:,end);
            else
                obj.LRKernel = lLRKernel(:,end);
                obj.HRKernel = lHRKernel(:,end);
            end
            obj.rf_indices_lr(end+1) = {IndexLR};
            obj.rf_indices_hr(end+1) = {IndexHR};
            obj.rf_models(end+1) = {'ADC'};
        end
        
        function obj = setData(obj,data)
            % Sets the target data for deconvolution. Please use this function
            % and do not assign the property directly.
            obj.data = data;
            obj.data_size = size(data.img);
            obj.LRKernel = [];
            obj.HRKernel = [];
            obj.n_anisotropic = 0;
            obj.n_isotropic = 0;
        end
        
        function obj = ReshapeDataToOriginal(obj)
            % After PerformDeconv is called, the data might be in vectorial
            % format. If you need it in the original shape, this allows for it
            obj.data = reshape(obj.data,obj.data_size);
        end
        
        function obj = AddAnisotropicRF_DKI(obj,EigenValues,MeanKurtosis)
            % Add an anisotropic RF based on the DTI/DKI model. These can be
            % added only before any isotropic RF. Eigenvalues is a vector
            % containing the three main eigenvalues of the diffusion tensor,
            % whereas meankurtosis is a scalar.
            if(isempty(obj.data))
                warning('Cannot initialize RFs if data has not been set. Bailing out.');
                return
            end
            if(obj.n_isotropic ~= 0)
                warning('Please add all anisotropic RFs before any isotropic component.');
                return
            end
            pass_data = obj.data;
            pass_data.bvecs(:,3) = -pass_data.bvecs(:,3);
            obj.n_anisotropic = obj.n_anisotropic+1;
            [~,lLRKernel,lHRKernel] = mDRLMT_MakeDKIKernel_multicomp(pass_data,obj.nDirections,EigenValues,MeanKurtosis,[],0);
            IndexHR = size(obj.HRKernel,2)+1:size(obj.HRKernel,2)+obj.nDirections;
            IndexLR = size(obj.LRKernel,2)+1:size(obj.LRKernel,2)+size(lLRKernel,2);
            if(~isempty(obj.LRKernel))
                obj.LRKernel(:,IndexLR) = lLRKernel;
                obj.HRKernel(:,IndexHR) = lHRKernel;
            else
                obj.LRKernel = lLRKernel;
                obj.HRKernel = lHRKernel;
            end
            obj.rf_indices_lr(end+1) = {IndexLR};
            obj.rf_indices_hr(end+1) = {IndexHR};
            obj.rf_models(end+1) = {'DKI'};
        end
        
        function obj = AddAnisotropicRF_NODDI(obj,noddi_parameters)
            % Add an anisotropic RF based on the NODDI model. These can be
            % added only before any isotropic RF.
            % noddi_parameters is a cell array containing elements (x), which should be specified as follows:
            % x(1) is the volume fraction of the intracellular space.
            % x(2) is the free diffusivity of the material inside and outside the cylinders.
            % x(3) is the concentration parameter of the Watson's distribution.
            % x(4) is the volume fraction of the isotropic compartment.
            % x(5) is the diffusivity of the isotropic compartment.
            % x(6) is the measurement at b=0.;
            
            if(isempty(obj.data))
                warning('Cannot initialize RFs if data has not been set. Bailing out.');
                return
            end
            if(obj.n_isotropic ~= 0)
                warning('Please add all anisotropic RFs before any isotropic component.');
                return
            end
            obj.n_anisotropic = obj.n_anisotropic+length(noddi_parameters);
            [~,lLRKernel,lHRKernel] = mDRLMT_MakeNODDIKernel_multicomp(obj.data,obj.nDirections,noddi_parameters,[],0);
            IndexHR = size(obj.HRKernel,2)+1:size(obj.HRKernel,2)+obj.nDirections;
            IndexLR = size(obj.LRKernel,2)+1:size(obj.LRKernel,2)+size(lLRKernel,2);
            if(~isempty(obj.LRKernel))
                obj.LRKernel(:,IndexLR) = lLRKernel;
                obj.HRKernel(:,IndexHR) = lHRKernel;
            else
                obj.LRKernel = lLRKernel;
                obj.HRKernel = lHRKernel;
            end
            obj.rf_indices_lr(end+1) = {IndexLR};
            obj.rf_indices_hr(end+1) = {IndexHR};
            obj.rf_models(end+1) = {'NODDI'};
        end
        
        function obj = AutomaticDRLDamping(obj)
            % Compute the automatic damping for the dRL method. Only needed if
            % the deconvolution method is dRL.
            obj.NN_L = get_drl_nn_heuristic(obj.LRKernel,max(obj.data.bvals),0.7e-3);
            obj.NN_H = get_drl_nn_heuristic(obj.HRKernel,max(obj.data.bvals),0.7e-3);
        end
        
        function output = PerformDeconv(obj,nof_workers,low_mem_mode)
            % Actually perform the deconvolution. All parameters must be set
            % before hand using the dedicated functions. Output is a structure
            % with fields FOD and FOD_norm. FOD_norm is a rescaled version of the FOD
            % meant to be compatible with ExploreDTI-based fiber tractography.
            
            if(nargin < 2)
                myCluster = parcluster('local');
                nof_workers = myCluster.NumWorkers;
            end
            
            if(nargin < 3)
                low_mem_mode = 0;
            end
            
            if(~obj.isInitialized())
                warning('Not all conditions are met to perform the deconvolution. Make sure you setted all needed parameters.');
                output = [];
                return
            end
            
            siz = obj.data_size;
            if(~ismatrix(obj.data.img))
                obj.data.img = reshape(obj.data.img,siz(1)*siz(2)*siz(3),siz(4)); % st(133)* ~isnan(FA)
                obj.data.img = permute(obj.data.img,[2 1]);
            end
            [st,sm] = size(obj.data.img);
            
            NC = obj.n_isotropic;
            ANC = obj.n_anisotropic;
            
            nreconstruction_vertices = size(obj.HRKernel,2)-NC;
            nreconstruction_vertices_lr = size(obj.LRKernel,2)-NC;
            
            %             fprintf('Determining NN: %.3f for LR and %.3f for HR %s',obj.NN_L, obj.NN_H, newline);
            
            shell_weight = zeros(size(obj.data.bvals));
            all_shells = (unique(obj.data.bvals));
            for ij=1:length(all_shells)
                this_shell = abs(obj.data.bvals-all_shells(ij)) < 100;
                shell_weight(this_shell) = ij;
            end
            
            shell_weight(shell_weight < length(all_shells)) = obj.inner_shells_weight;
            shell_weight(shell_weight == length(all_shells)) = 1;
            
            weighted_LRKernel = obj.LRKernel;
            weighted_HRKernel = obj.HRKernel;
            
            % Weighted kernels for the lower shells
            for ij=1:size(obj.LRKernel,2)
                weighted_LRKernel(:,ij) = obj.LRKernel(:,ij).*shell_weight;
            end
            
            for ij=1:size(obj.HRKernel,2)
                weighted_HRKernel(:,ij) = obj.HRKernel(:,ij).*shell_weight;
            end
            
            fractions = zeros([sm NC+ANC]);
            if(low_mem_mode == 0)
                WM_fod = zeros([sm nreconstruction_vertices],'single');
            else
                WM_fod = sparse([sm nreconstruction_vertices],'single');
            end
            RSS = zeros([sm 1]);
            S0 = zeros([sm 1]);
            
            N = sm;
            op_e2 = optimset('TolX',1e-2);
            
            [~,DeconvMethodCode] = SphericalDeconvolution.isSupportedMethod(obj.deconv_method); % -1 = failure; 1 = LSQNONNEG; 2 = DW_RegularizedDeconv; 3 = dRL
            if(DeconvMethodCode == -1)
                warning('Unsupported deconvolution method.');
                return;
            end
            
            TheTol = 1e-3;
            tic
            
            if(NC > 0)
                parfor (x=1:N,nof_workers)
                    %     if(mod(x,progressStepSize) == 0)
                    %         ppm.increment();
                    %     end
                    if(obj.data.mask(x) < 1)
                        continue
                    end

                    Stot = double(obj.data.img(:,x));

                    % This normalization assumes there is some b=0s/mm2 data. It is not
                    % essential to do it as long as fractions is normalized at the end
                    NormFactor = mean(Stot(obj.data.bvals<100));
                    S0(x) = NormFactor;
                    Stot = Stot/NormFactor;

                    Stot = Stot.*shell_weight; % Weight the lower shells

                    piso = zeros(NC,1);
                    p_old = Inf; % store the fractions at previous iter

                    % The following loop will iterate WM-FOD estimation (with mDRL) and
                    % fractions guess (with LSQNONNEG)
                    for iter=1:50 % 50 = max number of iterations but usually exits way before

                        DS = max(Stot-weighted_LRKernel(:,end-NC+1:end)*piso,0); % Subtract GM and CSF contributions

                        if(DeconvMethodCode == 4)
                            %                         fODFC = mat_dRL(DS, weighted_LRKernel(:,1:end-NC), 200, obj.NN_L, 8);
                            fODFC = ADT_deconv_RLdamp_1D_noEP(DS, weighted_LRKernel(:,1:end-NC),200, obj.NN_H);
                        elseif(DeconvMethodCode == 3)
                            %                         fODFC = mat_RL(DS, weighted_LRKernel(:,1:end-NC), 200);
                            fODF = RichardsonLucy(DS, weighted_LRKernel(:,1:end-NC), 200);
                        elseif(DeconvMethodCode == 2)
                            fODFC = DW_RegularizedDeconv(weighted_LRKernel(:,1:end-NC),DS,op_e2,obj.L2LSQ_reg);
                        elseif(DeconvMethodCode == 1)
                            fODFC = lsqnonneg(weighted_LRKernel(:,1:end-NC),DS,op_e2);
                        end
                        % This line is quite tricky. It comes out of some trial and error
                        % but actually has the intent of eliminating 1) small contributions
                        % 2) flat - spread fods enforcing sparsity.
                        % Without this line the code DOESN'T work. (WM is over-estimated).
                        fODFC(fODFC < median(fODFC)) = 0;

                        % Build a dictionary to fit the complete signal (Stot)
                        Y = [obj.LRKernel(:,1:end-NC)*fODFC obj.LRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)

                        if(sum(Y(:,1)) > 0) % i.e. if the FOD is non-zero
                            Y(:,1) = Y(:,1)/max(Y(:,1)); % Normalize WM signal
                        end

                        p = lsqnonneg(Y,Stot./shell_weight,op_e2); % Compute the signal fractions
                        piso = p(end-NC+1:end);

                        % if nothing changed compared to previous iter, exit. (Tol may need to be
                        % adjusted)
                        if(sum(abs(p-p_old) < TheTol) == 3)
                            break
                        end
                        p_old = p;
                    end

                    % NEW FINAL STEP 05/02/2018
                    DS = max(Stot - weighted_HRKernel(:,end-NC+1:end)*piso,0); % Subtract GM and CSF contributions

                    if(DeconvMethodCode == 4)
                        %                     fODF = mat_dRL(DS, weighted_HRKernel(:,1:end-NC),200, obj.NN_H, 8);
                        fODF = ADT_deconv_RLdamp_1D_noEP(DS, weighted_HRKernel(:,1:end-NC),200, obj.NN_H);
                    elseif(DeconvMethodCode == 3)
                        %                     fODF = mat_RL(DS, weighted_HRKernel(:,1:end-NC), 200);
                        fODF = RichardsonLucy(DS, weighted_HRKernel(:,1:end-NC), 200);
                    elseif(DeconvMethodCode == 2)
                        fODF = DW_RegularizedDeconv(weighted_HRKernel(:,1:end-NC),DS,op_e2, obj.L2LSQ_reg);
                    elseif(DeconvMethodCode == 1)
                        fODF = lsqnonneg(weighted_HRKernel(:,1:end-NC),DS, op_e2);
                    end
                    fODFC = fODF;
                    fODFC(fODFC < median(fODFC)) = 0;
                    if(obj.n_anisotropic == 1)
                        Y = [obj.HRKernel(:,1:end-NC)*fODFC obj.HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                        if(sum(Y(:,1)) > 0) % i.e. if the FOD is non-zero
                            Y(:,1) = Y(:,1)/max(Y(:,1)); % Normalize WM signal
                        end
                    else
                        Y = [];
                        cindex = 1;
                        for kANC = 1:ANC
                            Y = [Y weighted_HRKernel(:,cindex:cindex+obj.nDirections-1)*fODFC(cindex:cindex+obj.nDirections-1)];
                            cindex = cindex+obj.nDirections;
                            %         if(sum(Y(:,1)) > 0) % i.e. if the FOD is non-zero
                            % changed for multi FOD
                            if(sum(Y(:,kANC)) > 0) % i.e. if the FOD is non-zero
                                Y(:,kANC) = Y(:,kANC)/max(Y(:,kANC)); % Normalize WM signal
                            end
                        end
                        Y = [Y weighted_HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                    end
                    p = lsqnonneg(Y,Stot./shell_weight,op_e2); % Compute the signal fractions
                    RSS(x) = sum((Stot-Y*p).^2);

                    fractions(x,:) = p;
                    WM_fod(x,:) = single(fODF);
                end
                toc
            else
                parfor (x=1:N,nof_workers)
                    %     if(mod(x,progressStepSize) == 0)
                    %         ppm.increment();
                    %     end
                    if(obj.data.mask(x) < 1)
                        continue
                    end

                    Stot = double(obj.data.img(:,x));

                    % This normalization assumes there is some b=0s/mm2 data. It is not
                    % essential to do it as long as fractions is normalized at the end
                    NormFactor = mean(Stot(obj.data.bvals<100));
                    S0(x) = NormFactor;
                    Stot = Stot/NormFactor;
                    
                    fractions(x,1) = 1;
                    FOD = ADT_deconv_RLdamp_1D_noEP(Stot, weighted_HRKernel,200, obj.NN_H);
                    WM_fod(x,:) = single(FOD);
                    RSS(x) = sum((Stot-weighted_HRKernel*FOD).^2);
                end
            end
            % Restructure data
            
            fsum = sum(fractions,2);
            for ij=1:size(fractions,2)
                fractions(:,ij) = fractions(:,ij) ./ (fsum+eps);
            end
            
            output.fractions = reshape(fractions,[siz(1:3),size(fractions,2)]);
            
            WM_fod_max = max(WM_fod,[],2);
            WM_fod_normalized = WM_fod;
            WM_fod_val = mean(WM_fod_max(fractions(:,1) > 0.7*max(WM_fod_max(:)))); % 20/12/2017
            for ij=1:size(WM_fod_normalized,2)
                WM_fod_normalized(:,ij) = WM_fod_normalized(:,ij) / WM_fod_val;% .* fractions(:,1); % 20/12/2017
            end
            
            output.RSS = reshape(RSS,siz(1:3));
            clear RSS;
            output.FOD = reshape(WM_fod,[siz(1:3),size(WM_fod,2)]);
            clear WM_fod;
            output.FOD_norm = reshape(WM_fod_normalized,[siz(1:3),size(WM_fod_normalized,2)]);
            clear WM_fod_normalized;
            
            obj.data.img = permute(obj.data.img,[2 1]);
            obj.data.img = reshape(obj.data.img,siz); % st(133)* ~isnan(FA)
        end
        
        function setDeconvMethod(obj,method)
            % Sets the deconvolution method. one between csd, mscsd, grl, mfod
            if(~SphericalDeconvolution.isSupportedMethod(method))
                warning('Unsupported deconvolution method.');
                return;
            end
            obj.deconv_method = method;
        end
        
        % Sets the L2 regularization value.
        function setL2REG(obj,reg_val)
            obj.L2LSQ_reg = reg_val;
        end
        
    end
    
    methods(Static)
        function methods = SupportedMethods()
            % List the supported deconvolution methods
            methods = {'LSQ','L2LSQ','RL','dRL'};
        end
        
        % Check whether a method is actually supported
        function [boolean,method_id] = isSupportedMethod(method)
            methods = SphericalDeconvolution.SupportedMethods();
            for method_id=1:length(methods)
                if(strcmpi(methods{method_id},method))
                    boolean = true;
                    return;
                end
            end
            boolean = false;
            method_id = -1;
        end
        
        function SaveOutputToNii(SpherDec,output,file_prefix)
            % Save the content of a deconvolution data structure to nifti. SpherDec
            % is an instance of this class, output is the structure returned
            % from PerformDeconv. file_prefix is the name without extension of
            % the desired outputs.
            lmax = 16;
            super_scheme = gen_scheme(SpherDec.nDirections,lmax); % the reconstruction scheme. Change 300 to any number
            sh = SH(lmax,super_scheme.vert);
            
            if(SpherDec.n_anisotropic == 1)
                % Single FOD case
                fod = sh.coef(output.FOD);
                
                %             FOD_max = max(output.FOD,[],4);
                %             FOD_scaled = output.FOD;
                %             FOD_val = mean(FOD_max(output.fractions(:,:,:,1) > 0.7*max(FOD_max(:)))); % 20/12/2017
                %             for ij=1:size(FOD_scaled,4)
                %                 FOD_scaled(:,:,:,ij) = FOD_scaled(:,:,:,ij) / FOD_val;% .* fractions(:,:,:,1); % 20/12/2017
                %             end
                mult = 1;
                if(SpherDec.NN_H ~= 0)
                    mult = 0.1/SpherDec.NN_H;
                end
                
                if(isfield(SpherDec.data,'hdr'))
                    DW_SaveVolumeLikeNii(fod,SpherDec.data,[file_prefix '_CSD_FOD'],0);
                    DW_SaveVolumeLikeNii(output.fractions,SpherDec.data,[file_prefix '_fractions'],0);
                    DW_SaveVolumeLikeNii(FOD_scaled,SpherDec.data,[file_prefix '_CSD_FOD_scaled'],0);
                else
                    data_struct.img = fod;
                    data_struct.VD = SpherDec.data.VD;
                    EDTI.WriteNifti(data_struct,[file_prefix '_CSD_FOD.nii']);
                    data_struct.img = output.fractions;
                    EDTI.WriteNifti(data_struct,[file_prefix '_fractions.nii']);
                    data_struct.img = sh.coef(output.FOD*mult);
                    EDTI.WriteNifti(data_struct,[file_prefix '_CSD_FOD_scaled.nii']);
                end
            else
                % multi-FOD case
                for fod_id=1:SpherDec.n_anisotropic
                    fod = sh.coef(output.FOD(:,:,:,1+(fod_id-1)*SpherDec.nDirections:fod_id*SpherDec.nDirections));
                    
                    if(isfield(SpherDec.data,'hdr'))
                        DW_SaveVolumeLikeNii(fod,SpherDec.data,[file_prefix '_CSD_FOD_' num2str(fod_id)],0);
                    else
                        data_struct.img = fod;
                        data_struct.VD = SpherDec.data.VD;
                        EDTI.WriteNifti(data_struct,[file_prefix '_CSD_FOD_' num2str(fod_id) '.nii']);
                    end
                end
                if(isfield(SpherDec.data,'hdr'))
                    DW_SaveVolumeLikeNii(output.fractions,SpherDec.data,[file_prefix '_fractions'],0);
                else
                    data_struct.img = output.fractions;
                    EDTI.WriteNifti(data_struct,[file_prefix '_fractions.nii']);
                end
            end
        end
        
        function mrt_data = LoadNiiBvalBvec(varargin)
            % Load diffusion data in the MRIToolkit format. Input:
            % nii_file: the .nii(gz) data file
            % bval: the .bval file
            % bvec: the .bvec file
            % mask: the associated mask (optional)
            
            if( isempty(varargin))
                help('SphericalDeconvolution.LoadNiiBvalBvec');
                return;
            end
            
            coptions = varargin;
            file_in = GiveValueForName(coptions,'nii_file');
            if(isempty(file_in))
                error('Need to specify the input .nii file');
            end
            coptions = varargin;
            bval = GiveValueForName(coptions,'bval');
            if(isempty(bval))
                error('Need to specify the input .bval file');
            end
            coptions = varargin;
            bvec = GiveValueForName(coptions,'bvec');
            if(isempty(bvec))
                error('Need to specify the input .bvec file');
            end
            mask = GiveValueForName(coptions,'mask');
            if(isempty(mask))
                mask = '';
            end
            
            mrt_data = DW_LoadData(file_in,bvec,bval,mask);
            mrt_data.img = single(mrt_data.img);
            mrt_data.img = permute(mrt_data.img,[2 1 3 4]);
            mrt_data.img = flip(mrt_data.img,1);
            mrt_data.img = flip(mrt_data.img,2);
        end
        
        function [EigVal,IsoK] = EstimatedAverageEigval_IsoK(data)
            % Estimate the average eigenvalues and isotropic kurtosis in a
            % dataset to initialize the RF
            [sx,sy,sz,st] = size(data.img);
            signal_stack = reshape(data.img,sx*sy*sz,st);
            
            [G,WG,GoodIndexes] = DW_BuildDTMat(data,unique(data.bvals),1);
            Gt = G;
            G = [Gt(:,1:6) 1/6*(data.bvals*1e-3).^2 Gt(:,end)];
            
            msk = sum(signal_stack==0,2) == 0;
            p = zeros(length(signal_stack),8);
            K = zeros(length(signal_stack),1);
            autovals = zeros(length(signal_stack),3);
            parfor ij=1:length(signal_stack)
                if(msk(ij) < 1)
                    continue
                end
                S = signal_stack(ij,:)';
                if(sum(S==0) > 0)
                    continue
                end
                p1 = E_DTI_WLLS_WW(S,G)
                p(ij,:) = p1;%G\log(S);
                if(sum(~isfinite(p(ij,:))) > 0)
                    continue
                end
                autoval = eig(D2Dtensor(p(ij,:)));
                autovals(ij,:) = sort(autoval,'descend');
                K(ij) = 1e-6*p1(end-1)/mean(autoval.^2);
            end
            
            p = reshape(p,sx,sy,sz,8);
            
            [~,m_FA,m_DEC,~,~,m_lambdas,eigenvectors] = DW_ComputeTensorMetrics(p,1:6);
            U = m_FA(:) > 0.7 & data.mask(:) > 0;
            
            EigVal = mean(autovals(U > 0,:));
            EigVal(2:3) = mean(EigVal(2:3));
            IsoK = mean(K(U(:) > 0 & K(:) > 0 & K(:) < 4));
        end
        
        function TerminateTractsWithFraction(varargin)
            % This function filters a .MAT fiber tractography result using the
            % fractions estimated from the GRL method, to stop fiber
            % tractography at the GM-WM interface, or GM-CSF. This code assumes GRL
            % has been run with 3 classes (WM, GM, CSF). Input arguments
            % are:
            % mat_file: the reference ExploreDTI-like .MAT file
            % tract_file: the tractography file in ExploreDTI-like .MAT format
            % out_file: the desired output (.MAT)
            % fraction_file: 'The tissue class fractions (.nii)'
            % mask_mode: one between 'wm' (WM-GM interface) 'gm' (GM-CSF
            % interface) 'gm_only' (only the GM part of a tract)
            
            if( isempty(varargin))
                help('SphericalDeconvolution.TerminateTractsWithFraction');
                return;
            end
            
            coptions = varargin;
            mat_file = GiveValueForName(coptions,'mat_file');
            if(isempty(mat_file))
                error('Need to specify the reference .MAT file');
            end
            tract_file = GiveValueForName(coptions,'tract_file');
            if(isempty(tract_file))
                error('Need to specify the tracts .MAT file');
            end
            out_file = GiveValueForName(coptions,'out_file');
            if(isempty(out_file))
                error('Need to specify the output .MAT file');
            end
            fraction_file = GiveValueForName(coptions,'fraction_file');
            if(isempty(fraction_file))
                error('Need to specify the fractions .nii file');
            end
            mask_mode = GiveValueForName(coptions,'mask_mode');
            if(isempty(mask_mode))
                error('Need to specify the mask_mode');
            end
            
            load(tract_file);
            load(mat_file,'FA','VDims');
            
            [sx,sy,sz] = size(FA);
            
            intersect_mask = EDTI.LoadNifti(fraction_file);
            intersect_mask = intersect_mask.img;
            [~,intersect_mask] = max(intersect_mask,[],4);
            if(strcmp(mask_mode,'wm'))
                intersect_mask = intersect_mask == 1 & ~isnan(FA);
            elseif(strcmp(mask_mode,'gm'))
                intersect_mask = intersect_mask ~= 3 & intersect_mask > 0 & ~isnan(FA);
            elseif(strcmp(mask_mode,'gm_only'))
                intersect_mask = intersect_mask == 2 & ~isnan(FA);
            end
            
            new_tracts = {};
            new_tracts_L = {};
            new_tracts_FA = {};
            new_tracts_MD = {};
            new_tracts_FE = {};
            new_tracts_GEO = {};
            new_tracts_Lambdas = {};
            new_tracts_Angle = {};
            new_tracts_FOD = {};
            
            extra_tracts = 0;
            
            for tid=1:length(Tracts)
                tract = Tracts{tid};
                
                points2keep = false(size(tract,1),1);
                for l=1:size(tract,1)
                    point = round(tract(l,:)./VDims);
                    if(intersect_mask(point(1),point(2),point(3)))
                        % This is a good point
                        points2keep(l) = true;
                    end
                end
                
                if(sum(points2keep) > 0)
                    LabeledVector = LabelVector(points2keep);
                    for ij=2:max(LabeledVector)
                        points2keep_l = LabeledVector == ij;
                        if(sum(points2keep_l > 0) < 2)
                            continue
                        end
                        
                        TP = tract(points2keep_l,:);
                        TL = 0;
                        for tpid=2:size(TP,1)
                            TL = TL + norm(TP(tpid,:)-TP(tpid-1,:));
                        end
                        if(TL < Parameters.Length_range(1))
                            continue
                        end
                        
                        new_tracts{end+1} = TP;
                        new_tracts_L{end+1} = size(new_tracts{end},1);
                        new_tracts_FA{end+1} = TractFA{tid}(points2keep_l);
                        new_tracts_MD{end+1} = TractMD{tid}(points2keep_l);
                        new_tracts_GEO{end+1} = TractGEO{tid}(points2keep_l,:);
                        new_tracts_FE{end+1} = TractFE{tid}(points2keep_l,:);
                        new_tracts_Lambdas{end+1} = TractLambdas{tid}(points2keep_l,:);
                        new_tracts_FOD{end+1} = TractsFOD{tid}(points2keep_l);
                        new_tracts_Angle{end+1} = TractAng{tid}(points2keep_l);
                        extra_tracts = extra_tracts+1;
                        
                    end
                    points2keep(LabeledVector > 1) = false;
                end
                
                if(sum(points2keep) < 2)
                    points2keep = false(size(points2keep));
                end
                tract = tract(points2keep,:);
                Tracts{tid} = tract;
                TractL{tid} = size(tract,1);
                TractFA{tid} = TractFA{tid}(points2keep);
                TractMD{tid} = TractMD{tid}(points2keep);
                TractGEO{tid} = TractGEO{tid}(points2keep,:);
                TractFE{tid} = TractFE{tid}(points2keep,:);
                TractLambdas{tid} = TractLambdas{tid}(points2keep,:);
                TractsFOD{tid} = TractsFOD{tid}(points2keep);
                TractAng{tid} = TractAng{tid}(points2keep);
                
            end
            
            good_tracts = true(size(Tracts));
            for tract_id=1:length(Tracts)
                if(TractL{tract_id} == 0)
                    good_tracts(tract_id) = false;
                end
            end
            
            Tracts = Tracts(good_tracts);
            TractsFOD = TractsFOD(good_tracts);
            TractL = TractL(good_tracts);
            TractFA = TractFA(good_tracts);
            TractFE = TractFE(good_tracts);
            TractAng = TractAng(good_tracts);
            TractGEO = TractGEO(good_tracts);
            TractLambdas = TractLambdas(good_tracts);
            TractMD = TractMD(good_tracts);
            
            
            Tracts(end+1:end+extra_tracts) = new_tracts;
            TractL(end+1:end+extra_tracts) = new_tracts_L;
            TractFA(end+1:end+extra_tracts) = new_tracts_FA;
            TractMD(end+1:end+extra_tracts) = new_tracts_MD;
            TractFE(end+1:end+extra_tracts) = new_tracts_FE;
            TractGEO(end+1:end+extra_tracts) = new_tracts_GEO;
            TractLambdas(end+1:end+extra_tracts) = new_tracts_Lambdas;
            TractAng(end+1:end+extra_tracts) = new_tracts_Angle;
            TractsFOD(end+1:end+extra_tracts) = new_tracts_FOD;
            FList = (1:length(Tracts))';
            
            % clear new_tracts new_tracts_L new_tracts_FA new_tracts_MD new_tracts_FE new_tracts_GEO
            % clear new_tracts_Lambdas new_tracts_Angle new_tracts_FOD points2keep tracts
            
            %             disp(['Gated ' num2str(sum(good_tracts == true)) ' out of ' num2str(length(good_tracts)) ' ' ...
            %                 sprintf('%.2f',100*single(sum(good_tracts == true))/single(length(good_tracts))) '%']);
            
            save(out_file,'Tracts','TractsFOD','TractL','TractFA','TractFE','TractFE',...
                'TractAng','TractGEO','TractLambdas','TractMD','FList','TractMask','VDims','-v7.3');
        end
        
        function [init_lambdas,init_K] = Eigenval_IsoK_WM_FromData(data,FA,mask)
            % Estimate the tensor eigenvalues and isotropic kurtosis from the data, in a mask where FA >= 0.7
            data.img = single(data.img);
            [sx,sy,sz,st] = size(data.img);
            signal_stack = reshape(data.img,sx*sy*sz,st);
            
            [G,WG,GoodIndexes] = DW_BuildDTMat(data,unique(data.bvals),1);
            Gt = G;
            G = [Gt(:,1:6) 1/6*(data.bvals*1e-3).^2 Gt(:,end)]; % Extend for isotropic kurtosis
            msk = sum(signal_stack==0,2) == 0;
            p = zeros(length(signal_stack),8);
            K = zeros(length(signal_stack),1);
            autovals = zeros(length(signal_stack),3);
            parfor ij=1:length(signal_stack)
                if(msk(ij) < 1)
                    continue
                end
                S = signal_stack(ij,:)';
                if(sum(S==0) > 0)
                    continue
                end
                p1 = E_DTI_WLLS_WW(S,G); % This uses the weighted least squares as implemented in ExploreDTI
                p(ij,:) = p1;%G\log(S);
                if(sum(~isfinite(p(ij,:))) > 0)
                    continue
                end
                autoval = eig(D2Dtensor(p(ij,:)));
                autovals(ij,:) = sort(autoval,'descend');
                K(ij) = 1e-6*p1(end-1)/mean(autoval.^2);
            end
            
            p = reshape(p,sx,sy,sz,8);
            
            U = FA(:) > 0.7 & mask(:) > 0;
            
            init_lambdas = real(mean(autovals(U > 0,:)));
            init_lambdas(2:3) = mean(init_lambdas(2:3));
            init_K = real(mean(K(U(:) > 0 & K(:) > 0 & K(:) < 4)));
            
            disp(['Calibrated lambdas are:' num2str(init_lambdas)]);
            disp(['Calibrated K is:' num2str(init_K)]);
            
            data.img = reshape(data.img,sx,sy,sz,st);
            
        end
        
        function sh_matrix = SHFitMatrix(varargin)
            % constructs the Spherical Harmonics fit matrices used in
            % ExploreDTI, Dipy, mrtrix3.
            % Input arguments:
            % bvecs: the gradient orientations of the data on the unit
            % sphere (file or vector)
            % bvals: the corresponding diffusion weightings (file or
            % matrix)
            % bmat: exclusive with bvals/bvecs. Use the b-matrix from
            % ExploreDTI instead (file or matrix)
            % basis: one in "edti" (default), "dipy" (descoteaux et al.),
            % "mrtrix" (tournier et al.)
            % bvalue: the b-value to fit. If not specified, defaults to the
            % maximum shell
            % lmax: the maximum spherical harmonics order. default 8 or
            % what is allowed by the protocol
            bvals_file = GiveValueForName(varargin,'bvals');
            bvecs_file = GiveValueForName(varargin,'bvecs');
            txt_file = GiveValueForName(varargin,'bmat');
            basis_type = GiveValueForName(varargin,'basis');
            bvalue = GiveValueForName(varargin,'bvalue');
            lmax = GiveValueForName(varargin,'lmax');

            if(isempty(txt_file) && (isempty(bvals_file) || isempty(bvecs_file)))
                error('Missing mandatory argument bvals/bvecs or txt');
            end
            
            if(isempty(basis_type))
                basis_type = 'edti';
            end
            
            if(~isempty(txt_file))
               if(ischar(txt_file))
                   bmat = load(txt_file);
               else
                   bmat = txt_file;
               end
               [bvals,bvecs] = EDTI.bval_bvec_from_b_Matrix(bmat); 
            else
                if(ischar(bvals_file))
                    bvals = load(bvals_file);
                    bvals = bvals';
                else
                    bvals = bvals_file;
                end
                if(ischar(bvecs_file))
                    bvecs = load(bvecs_file);
                    bvecs = bvecs';
                else
                    bvecs = bvecs_file;
                end
            end
            
            if(isempty(bvalue))
                bvalue = max(round(bvals));
            end
            
            IX = abs(bvals-bvalue) < 0.1*bvalue;
%             bvals = bvals(IX);
            bvecs = bvecs(IX,:);
                        
            if(isempty(lmax))
                lmax = min(2.*(floor((sqrt(1+8.*size(bvecs,1))-3)./4)),32);
            end
            LmaxCoeffs = (lmax+1).*(lmax+2)./2;
            
            angle_rep = c2s(bvecs);
            sh_matrix = zeros(size(bvecs,1),LmaxCoeffs);
            if(strcmp(basis_type,'dipy'))
                bvecs_dipy = bvecs(:,[2 1 3]);
                bvecs_dipy(:,1) = -bvecs_dipy(:,1);
                angle_rep_dipy = c2s(bvecs_dipy);
                for l_order=0:2:lmax
                    order_legendre = legendre(l_order,cos(angle_rep_dipy(:,1)),'sch')';
                    for m_phase=-l_order:l_order
                        Ylm = (-1)^m_phase*sqrt((2*l_order+1)/(4*pi));
                        Ylm = Ylm.*order_legendre(:,1+abs(m_phase)).*exp(1i*m_phase*angle_rep_dipy(:,2));
                        index = (l_order^2+l_order+2)/2+m_phase;
                        switch(sign(m_phase))
                            case 0
                                sh_matrix(:,index) = Ylm;
                            case 1
                                sh_matrix(:,index) = imag(Ylm);
                            case -1
                                sh_matrix(:,index) = real(Ylm);
                        end
                    end
                end
            elseif(strcmp(basis_type,'edti'))
                for l_order=0:2:lmax
                    index = l_order*(l_order+1)/2 + 1;
                    order_legendre = legendre(l_order,cos(angle_rep(:,1)),'sch')';
                    for m_phase=0:l_order
                        Ylm = (-1)^m_phase*sqrt((2*l_order+1)/(4*pi));
                        Ylm = Ylm.*order_legendre(:,1+abs(m_phase)).*exp(1i*m_phase*angle_rep(:,2));
                        if(m_phase == 0)
                            sh_matrix(:,index) = Ylm;
                        else
                            sh_matrix(:,index-m_phase) = imag(Ylm);
                            sh_matrix(:,index+m_phase) = real(Ylm);
                        end
                    end
                end
            elseif(strcmp(basis_type,'mrtrix'))
%                 bvecs_mrtrix = bvecs;
                bvecs_mrtrix = bvecs(:,[2 1 3]);
                bvecs_mrtrix(:,1) = -bvecs_mrtrix(:,1);
                angle_rep_mrtrix = c2s(bvecs_mrtrix);
                for l_order=0:2:lmax
                    index = l_order*(l_order+1)/2 + 1;
                    order_legendre = legendre(l_order,cos(angle_rep_mrtrix(:,1)),'sch')';
                    for m_phase=0:l_order
                        Ylm = (-1)^m_phase*sqrt((2*l_order+1)/(4*pi));
                        Ylm = Ylm.*order_legendre(:,1+abs(m_phase)).*exp(1i*m_phase*angle_rep_mrtrix(:,2));
                        if(m_phase == 0)
                            sh_matrix(:,index) = Ylm;
                        else
                            sh_matrix(:,index-m_phase) = imag(Ylm);
                            sh_matrix(:,index+m_phase) = real(Ylm);
                        end
                    end
                end
            else
                error('Unknown basis type');
            end

        end
        
        function SHFittedData = dMRI_2_SH(varargin)
            % Converts dMRI data to the SH representation of ExploreDTI,
            % Dipy or mrtrix3.
            % Input arguments:
            % data: the input data structure. 
            % basis: the desired SH basis, one in
            % "edti" (default),"dipy","mrtrix"
            % bvalue: which b-value to fit. Defaults to the maximum
            % lmax: the order of the spherical harmonics
            data = GiveValueForName(varargin,'data');
            basis = GiveValueForName(varargin,'basis');
            bvalue = GiveValueForName(varargin,'bvalue');
            lmax = GiveValueForName(varargin,'lmax');
            if(isempty(data))
                error('Missing mandatory input data');
            end
            if(isempty(basis))
                basis = 'edti';
            end
            if(isempty(bvalue))
                bvalue = max(round(data.bvals));
            end      
            
            IX = abs(data.bvals-bvalue)<0.1*bvalue;
            data.bvals = data.bvals(IX);
            data.bvecs = data.bvecs(IX,:);
            if(~ismatrix(data.img))
                data.img = data.img(:,:,:,IX);
            else
                data.img = data.img(IX);
            end
            
            [sx,sy,sz,st] = size(data.img);
            if(~ismatrix(data.img))
                data.img = reshape(data.img,sx*sy*sz,st);
            else
                data.img = data.img';
            end
            SHB = SphericalDeconvolution.SHFitMatrix('basis',basis,'bvals',data.bvals,...
                'bvecs',data.bvecs,'lmax',lmax);
            SHFittedData = SHB\data.img';
            if(sx ~= 1 && sy ~= 1)
                SHFittedData = reshape(SHFittedData',sx,sy,sz,size(SHFittedData,1));
            end
        end
        
        function data = SH_2_dMRI(varargin)
            % Converts ExploreDTI, Dipy or mrtrix3 SH coefficients to
            % another program. Input arguments:
            % shcoeffs: the input representation
            % basis: the SH basis type, one in
            % "edti" (default),"dipy","mrtrix"
            % bvecs: the gradient orientations (variable or file)
            % bvals: the corresponding diffusion weightings (variable or
            % file)
            % bmat: alternatively, the b-matrix (variable or file)    
            % bvalue: the shell b-value. Defaults to the max
            % lmax: the order of the spherical harmonics
            bvals = GiveValueForName(varargin,'bvals');
            bvecs = GiveValueForName(varargin,'bvecs');
            basis = GiveValueForName(varargin,'basis');
            bvalue = GiveValueForName(varargin,'bvalue');
            bmat = GiveValueForName(varargin,'bmat');
            lmax = GiveValueForName(varargin,'lmax');
            shcoeffs = GiveValueForName(varargin,'shcoeffs');
            if(isempty(shcoeffs))
                error('Missing mandatory input shcoeffs');
            end
            
            if(isempty(basis))
                basis = 'edti';
            end
            
            if(isempty(bvals))
               [bvals,bvecs] = EDTI.bval_bvec_from_b_Matrix(bmat); 
            end
            
            if(ischar(bvals))
                bvals = load(bvals);
                bvecs = load(bvecs);
                bvals = bvals';
                bvecs = bvecs';
            end                
            
            if(isempty(bvalue))
                bvalue = max(round(bvals));
            end      
            
            IX = abs(bvals-bvalue)<0.1*bvalue;
            bvals = bvals(IX);
            bvecs = bvecs(IX,:);
            
            SHB = SphericalDeconvolution.SHFitMatrix('basis',basis,'bvals',bvals,...
                'bvecs',bvecs,'bmat',bmat,'lmax',lmax);

            [sx,sy,sz,st] = size(shcoeffs);
            if(~ismatrix(shcoeffs))
                shcoeffs = reshape(shcoeffs,sx*sy*sz,st)';
            end
            data.img = SHB*shcoeffs;           

            data.bvals = bvals;
            data.bvecs = bvecs;
            
            if(sx ~= 1 && sy ~= 1)
                data.img = data.img';
                data.img = reshape(data.img,sx,sy,sz,size(bvals,1));
            end
        end
        
        function SHCoeffs = SH_2_SH(varargin)
            % Converts ExploreDTI, Dipy or mrtrix3 SH coefficients to
            % another program. Input arguments:
            % shcoeffs: the input representation
            % basis_in: the input basis type, one in
            % "edti","dipy","mrtrix"
            % basis_out: the output basis type
            % bvecs: the gradient orientations (variable or file)
            % bvals: the corresponding diffusion weightings (variable or
            % file)
            % bmat: alternatively, the b-matrix (variable or file)
            % lmax: the order of the spherical harmonics
            % output: the output .nii file (optional) - only valid if
            % specifying the input as file
            bvals = GiveValueForName(varargin,'bvals');
            bvecs = GiveValueForName(varargin,'bvecs');
            basis_in = GiveValueForName(varargin,'basis_in');
            basis_out = GiveValueForName(varargin,'basis_out');
            bmat = GiveValueForName(varargin,'bmat');
            shcoeffs = GiveValueForName(varargin,'shcoeffs');
            lmax = GiveValueForName(varargin,'lmax');
            output = GiveValueForName(varargin,'output');
            if(isempty(shcoeffs))
                error('Missing mandatory input shcoeffs');
            end
            
            VD = [];
            if(ischar(shcoeffs))
                img = EDTI.LoadNifti(shcoeffs);
                shcoeffs = img.img;
                VD = img.VD;
                clear img;
            end
            
            if(isempty(bvals) && isempty(bmat))
               bvals = 3000*ones(300,1); 
               Q = gen_scheme(300,4);
               bvecs = Q.vert;
            end
            
            if(isempty(basis_in))
                error('Missing mandatory argument basis_in');
            end
            if(isempty(basis_out))
                error('Missing mandatory argument basis_out');
            end
            
            data = SphericalDeconvolution.SH_2_dMRI('bvals',bvals,...
                'bvecs',bvecs,'basis',basis_in,'bmat',bmat,'shcoeffs',shcoeffs,'lmax',lmax);
            SHCoeffs = SphericalDeconvolution.dMRI_2_SH('data',data,'basis',...
                basis_out,'lmax',lmax);
            
            if(~isempty(VD) && ~isempty(output))
               out.VD = VD;
               out.img = SHCoeffs;
               clear data SHCoeffs
               EDTI.WriteNifti(out,output); 
            end
        end
        
        function [fod_amp, vertices] = SH_FOD_2_Sphere(varargin)
            % Projects an FOD in Spherical Harmonics (in ExploreDTI basis) on the unit sphere
            % input arguments:
            % fod: the input FOD or the corresponding file
            % nvertices: the number of desired vertices
            % lmax: the spherical harmonics order
            fod_file = GiveValueForName(varargin,'fod');
            if(isempty(fod_file))
                error('Missing mandatory argument fod');
            end
            nvert = GiveValueForName(varargin,'nvertices');
            if(isempty(nvert))
                nvert = 720;
            end

            if(ischar(fod_file))
                DATA = EDTI.LoadNifti(fod_file);
            else
                DATA.img = fod_file;
                clear fod_file
            end
            lmax = GiveValueForName(varargin,'lmax');
            if(isempty(lmax))
                lmax = 16;
                switch(size(DATA.img,4))
                    case 15
                        lmax = 4;
                    case 28
                        lmax = 6;
                    case 45
                        lmax = 8;
                    case 163
                        lmax = 16;
                end
            end            
            
            SV = gen_scheme(nvert,4);
            vertices = SV.vert;
            
            mySH = SH(lmax,SV.vert);
            fod_amp = mySH.amp(DATA.img);
            
        end
        
        function [peak_dir,peak_amp,AFD,NuFO] = FODPeaksAndMetrics(varargin)
            % Determine the FOD peaks and their amplitude. Can also return
            % the AFD and NuFO metrics
            % input arguments:
            % fod: the input FOD or file reference
            % lmax: the SH order (default 8)
            % peak_threshold: minimum amplitude of a "genuine" peak
            % (default 0.1)
            % max_peaks: maximum allowed number of peaks (default 3)
            % parallel: allow parallel computing (default 1)
            fod = GiveValueForName(varargin,'fod');
            if(isempty(fod))
                error('Missing mandatory argument fod');
            end
            if(ischar(fod))
                C = EDTI.LoadNifti(fod);
                fod = C.img;
                clear C
            end
            lmax = GiveValueForName(varargin,'lmax');
            if(isempty(lmax))
                lmax = 8;
            end
            max_peaks = GiveValueForName(varargin,'max_peaks');
            if(isempty(max_peaks))
                max_peaks = 3;
            end
            peak_threshold = GiveValueForName(varargin,'peak_threshold');
            if(isempty(peak_threshold))
                peak_threshold = 0.1;
            end
                      
            SHPrecomp.init(lmax,300);
            c = parcluster('local');
            max_cores = c.NumWorkers;
            if(parallel == 0)
                max_cores = 1;
            end
            
            [peak_dir,peak_amp] = SHPrecomp.all_peaks(fod,peak_threshold,max_peaks,max_cores,lmax,300);
            
            if(nargout > 2)
               NuFO = zeros(size(peak_dir));
               AFD = zeros(size(peak_dir));
                           
                for x=1:size(AFD,1)
                    for y=1:size(AFD,2)
                        for z=1:size(AFD,3)
                            if(isempty(peak_dir{x,y,z}))
                                continue
                            end
                            NuFO(x,y,z) = size(peak_dir{x,y,z},2);
                            AFD(x,y,z) = peak_amp{x,y,z}(1);
                        end
                    end
                end
               
            end
            
        end
        
        function ConvertMatTractography2TCK(varargin)
            % Converts ExploreDTI tractography 2 the TCK format
            % input arguments:
            % mat_file: the input tractography mat file
            % output: the output tck file
            mat_file = GiveValueForName(varargin,'mat_file');
            if(isempty(mat_file))
                error('Missing mandatory argument mat_file');
            end
            output = GiveValueForName(varargin,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
           
            T = load(mat_file,'Tracts','TractMask','VDims');
            tracks.max_num_tracks = num2str(length(T.Tracts));
            tracks.unidirectional = 0;
            tracks.data = cell(1,length(T.Tracts));
            shift = round(size(T.TractMask).*T.VDims)/2;
            for ij=1:length(T.Tracts)
               TR = T.Tracts{ij}(:,[2 1 3]);
               TR(:,1) = -TR(:,1)+shift(2);
               TR(:,2) = -TR(:,2)+shift(1);
               TR(:,3) = TR(:,3) - shift(3);
               tracks.data(ij) = {TR}; 
            end
            write_mrtrix_tracks(tracks,output);
        end
     
        function ConvertTCKTractography2Mat(varargin)
            % Converts TCK tractography 2 the ExploreDTI format
            % input arguments:
            % tck_file: the input tractography mat file
            % nii_file: a reference .nii file in the same space of the tck
            % output: the output MAT file            
            tck_file = GiveValueForName(varargin,'tck_file');
            if(isempty(tck_file))
                error('Missing mandatory argument tck_file');
            end
            output = GiveValueForName(varargin,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
            nii_file = GiveValueForName(varargin,'nii_file');            
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            tracks = read_mrtrix_tracks(tck_file);
            tracks.max_num_tracks = length(tracks.data);
            Tracts = cell(1,tracks.max_num_tracks);
            TractL = cell(1,tracks.max_num_tracks);
            TractFA = cell(1,tracks.max_num_tracks);
            TractMD = cell(1,tracks.max_num_tracks);
            TractFE = cell(1,tracks.max_num_tracks);
            TractGEO = cell(1,tracks.max_num_tracks);
            TractAng = cell(1,tracks.max_num_tracks);
            TractLambdas = cell(1,tracks.max_num_tracks);
            FList = (1:length(Tracts))';
            
            ref = EDTI.LoadNifti(nii_file);
            VDims = ref.VD;
            TractMask = ones(size(ref.img));
            TractMask = permute(TractMask,[2 1 3]);
            shift = round(size(TractMask).*VDims)/2;
            for tid=1:length(tracks.data)
               T = tracks.data{tid}; 
               T(:,3) = T(:,3) + shift(3);
               T(:,2) = -T(:,2) + shift(2);
               T(:,1) = -T(:,1) + shift(1);
               T = T(:,[2 1 3]);
               
               Tracts(tid) = {T};
               TractFA(tid) = {ones(size(T,1),1)};
               TractMD(tid) = {ones(size(T,1),1)};
               TractAng(tid) = {ones(size(T,1),1)};
               TractFE(tid) = {ones(size(T,1),3)};
               for pid=2:size(T,1)
                  D = abs(T(pid,:)-T(pid-1,:));
                  D = D/norm(D);
                  TractFE{tid}(pid,:) = (D);
               end
               TractGEO(tid) = {ones(size(T,1),1)};
               TractLambdas(tid) = {ones(size(T,1),3)};
               TractL(tid) = {size(T,1)};
            end
            
            save(output,'Tracts','TractL','TractFA','TractFE','TractFE',...
                'TractAng','TractGEO','TractLambdas','TractMD','FList','TractMask','VDims','-v7.3');
            
        end
        
        function [mx,my,mz] = ConvertMatTractographt2VTK(varargin)
            % Converts MAT tractography 2 the VTK format
            % The output mx,my,mz are the shifts from the EDTI space to the
            % Slicer space
            % input arguments:
            % mat_file: the input tractography mat file
            % output: the output vtk file            
            mat_file = GiveValueForName(varargin,'mat_file');
            if(isempty(mat_file))
                error('Missing mandatory argument mat_file');
            end
            output = GiveValueForName(varargin,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
            [mx,my,mz] = E_DTI_Convert_tracts_mat_2_vtk_lines(mat_file, output);
        end
        
        function VTKTractography2Mat(varargin)
            % Converts VTK tractography 2 the ExploreDTI format
            % input arguments:
            % vtk_file: the input tractography mat file
            % nii_file: a reference .nii file in the same space of the tck
            % output: the output MAT file
            
            vtk_file = GiveValueForName(varargin,'vtk_file');
            if(isempty(vtk_file))
                error('Missing mandatory argument vtk_file');
            end
            output = GiveValueForName(varargin,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
            nii_file = GiveValueForName(varargin,'nii_file');            
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end
            
            REF = EDTI.LoadNifti(nii_file);
            
            try
                [Points,Lines] = VTKImport_test(fullfile(vtk_files(vid).folder,vtk_files(vid).name));
            catch 
                disp('Error trying to read the VTK file');
                return
            end
            Points(:,1) = Points(:,1) + mx;
            Points(:,2) = Points(:,2) + my;
            Points(:,3) = Points(:,3) + mz;
            Points(:,1) = size(REF.img,2) - Points(:,1); 
            Points(:,2) = size(REF.img,1) - Points(:,2); 
            Points = Points(:,[2 1 3]); 
            Tracts.VDims = REF.VD;
            Tracts.Tracts = cell(1,length(Lines));
            Tracts.TractAng = cell(1,length(Lines));
            Tracts.TractFA = cell(1,length(Lines));
            Tracts.TractFE = cell(1,length(Lines));
            Tracts.TractGEO = cell(1,length(Lines));
            Tracts.TractL = cell(1,length(Lines));
            Tracts.TractLambdas = cell(1,length(Lines));
            Tracts.TractMD = cell(1,length(Lines));
            Tracts.TractMask = ones(size(REF.img));
            Tracts.TractMask = permute(Tracts.TractMask,[2 1 3]);
            Tracts.FList = (1:length(Tracts.Tracts))';
            for tid=1:length(Tracts.Tracts)
                PL = Points(Lines{tid}+1,:);
                Tracts.Tracts(tid) = {PL};
                cTractFA = zeros([size(PL,1),1]);
                cTractFE = zeros([size(PL,1),3]);
                cTractGEO = zeros([size(PL,1),1]);
                cTractAng = zeros([size(PL,1),1]);
                cTractL = size(PL,1);
                cTractLambdas = zeros([size(PL,1),3]);
                cTractMD = zeros([size(PL,1),1]);
                for point=1:size(PL,1)
                   P = PL(point,:)./REF.VD; 
                   PC = floor(P)+1;
                   cTractFA(point) = DIFF_PROP.FA(PC(1),PC(2),PC(3))/sqrt(3);
                   cTractFE(point,:) = abs(squeeze(DIFF_PROP.FE(PC(1),PC(2),PC(3),:)));
                   cTractLambdas(point,:) = squeeze(DIFF_PROP.eigval(PC(1),PC(2),PC(3),:));
                   cTractMD(point) = mean(cTractLambdas(point,:));
                end
                Tracts.TractAng(tid) = {cTractAng};
                Tracts.TractFA(tid) = {cTractFA};
                Tracts.TractFE(tid) = {cTractFE};
                Tracts.TractGEO(tid) = {cTractGEO};
                Tracts.TractLambdas(tid) = {cTractLambdas};
                Tracts.TractMD(tid) = {cTractMD};       
                Tracts.TractL(tid) = {cTractL};
            end
            save(output,'-struct','Tracts');

        end
        
    end
end

% Private function to parse input parameters
function keys = ParseInputKeys(input)
keys = {};
if(mod(length(input),2) ~= 0)
    warning('Incorrect input pairs');
    return;
end
keys = cell(length(input)/2,2);
for key_id=1:2:length(input)
    keys(key_id,1) = {input{key_id}};
    keys(key_id,2) = {input{key_id+1}};
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

% Helper function
function LabeledVector = LabelVector(LabeledVector)
LabeledVector = uint32(LabeledVector);
CLabel = 1;
for ij=2:length(LabeledVector)
    if(LabeledVector(ij) > 0 && LabeledVector(ij-1) == 0)
        CLabel = CLabel + 1;
    end
    LabeledVector(ij) = LabeledVector(ij)*CLabel;
end
end

% From ExploreDTI: Weighted least squares
function X = E_DTI_WLLS_WW(S0,B)

warning off all


crit = 1;
iter_max = 10;
cn=0;
con = 10^-4;

Log_S0 = log(S0);

W = ones(length(S0),1);



X = ones(size(B,2),1);


while crit==1
    
    cn = cn+1;
    W_p = W;
    X_p = X;
    w = diag(W_p.^2);
    X = (B'*w*B)\(B'*w)*Log_S0;
    W = exp(B*X);
    
    if all(abs(X-X_p) <= con*max(abs(X),abs(X_p))) || cn==iter_max
        crit=0;
    end
    
end

if any(isnan(X)) || any(isinf(X))
    
    X = B\S0;
    
end


% disp(num2str(cn))

warning on all
end

% Private function to select a subset of b-values
function fitted_indexes = DW_GetFittedIndexes(data,fitted_bvalues)

fitted_indexes = [];
for j=1:length(fitted_bvalues)
    fitted_indexes = [fitted_indexes;ascolumn(find(data.bvals==fitted_bvalues(j)))];
end

end

% Generate a b-matrix from bval/bvec
function [G,WG,GoodIndexes] = DW_BuildDTMat(data,bvalues,fit_S0)
if(nargin < 2)
    bvalues_range = [min(data.bvals(data.bvals>=1)) max(data.bvals)];
end
GoodIndexes = DW_GetFittedIndexes(data,bvalues);
G = [data.bvecs(GoodIndexes,1).^2 data.bvecs(GoodIndexes,1).*data.bvecs(GoodIndexes,2) data.bvecs(GoodIndexes,1).*data.bvecs(GoodIndexes,3) ...
    data.bvecs(GoodIndexes,2).^2 data.bvecs(GoodIndexes,2).*data.bvecs(GoodIndexes,3) data.bvecs(GoodIndexes,3).^2];
%     G = [data.bvecs(GoodIndexes,1).^2 2*data.bvecs(GoodIndexes,1).*data.bvecs(GoodIndexes,2) 2*data.bvecs(GoodIndexes,1).*data.bvecs(GoodIndexes,3) ...
%         data.bvecs(GoodIndexes,2).^2 2*data.bvecs(GoodIndexes,2).*data.bvecs(GoodIndexes,3) data.bvecs(GoodIndexes,3).^2];
for j=1:6
    G(:,j) = -G(:,j).*data.bvals(GoodIndexes);
end
if(nargin > 2 && fit_S0 == 1)
    %         G = [G data.bvals(GoodIndexes)==0];
    G = [G ones(size(G,1),1)];
end
%     bvals = data.bvals(GoodIndexes);
%     bvecs = data.bvecs(GoodIndexes,:);
if(nargout > 1)
    WG = G;
    for h=1:6%size(WG,2)
        WG(:,h) = WG(:,h).*(data.bvals(GoodIndexes));
        WG(isnan(WG)) = 0;
        WG(~isfinite(WG)) = 0;
    end
end
end

% Private function assembling the DTI/DKI kernels
function [bmat,LRKernel,HRKernel] = mDRLMT_MakeDKIKernel_multicomp(data,nreconstruction_vertices,lambdas,K,isoDs,shell_data)

shells = single(unique(int16(round(data.bvals)/1)*1)); % automatic detection of number of shells.
% This shell splitting works only for b-value spaced more than 100. To be fixed for other datasets.
ndirections = zeros(length(shells),1);
for ij=1:length(shells)
    ndirections(ij) = sum(abs(data.bvals-shells(ij))<1);
end


if(shell_data > 0)
    bmat = cell(length(shells),1); % bmat = [bvecs bvals]
    Kernel = cell(length(shells),1);
    Kernel_LR = cell(length(shells),1);
else
    bmat{1} = [data.bvecs data.bvals];
    Kernel = cell(1,1);
    Kernel_LR = cell(1,1);
end

super_scheme = gen_scheme(nreconstruction_vertices,4); % the reconstruction scheme. Change 300 to any number
HRKernel = zeros(sum(ndirections),nreconstruction_vertices+length(isoDs));
[phi, theta] = cart2sph(super_scheme.vert(:,1),super_scheme.vert(:,2),super_scheme.vert(:,3)); % polar decomposition
lr_scheme = gen_scheme(min(length(data.bvals),90),4);
[phi_LR, theta_LR] = cart2sph(lr_scheme.vert(:,1),lr_scheme.vert(:,2),lr_scheme.vert(:,3));
LRKernel = zeros(sum(ndirections),size(lr_scheme.vert,1)+length(isoDs));

if(shell_data > 0)
    
    bvals = zeros(sum(ndirections),1);
    
    index = 1;
    
    for ij=1:length(shells)
        bmat{ij} = zeros(ndirections(ij),4);
        bmat{ij}(:,4) = shells(ij);
        bmat{ij}(:,1:3) = data.bvecs(index:index+ndirections(ij)-1,:);
        
        % Here the deconvolution dictionary is actually built
        Kernel{ij} = zeros(ndirections(ij),length(phi));
        for i=1:length(phi)
            anglesFi = [phi(i), theta(i)]*(180/pi); % in degrees
            Kernel{ij}(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K);
        end
        for l=1:length(isoDs)
            Kernel{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end
        
        Kernel_LR{ij} = zeros(ndirections(ij),length(phi_LR));
        for i=1:length(phi_LR)
            anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
            Kernel_LR{ij}(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K);
        end
        for l=1:length(isoDs)
            Kernel_LR{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end
        
        bvals(index:index+ndirections(ij)-1) = bmat{ij}(:,4);
        HRKernel(index:index+ndirections(ij)-1,:) = Kernel{ij};
        LRKernel(index:index+ndirections(ij)-1,:) = Kernel_LR{ij};
        
        % Just for the simulated signal - not needed in production version
        %    R = eye(3,3);
        %    D = diag(lambdas);
        %    S_WM = exp(-bmat{ij}(:,4).*diag(bmat{ij}(:,1:3)*R*D*R'*bmat{ij}(:,1:3)'));
        %    S_GM = exp(-bmat{ij}(:,4)*D_gm);
        %    S_CSF = exp(-bmat{ij}(:,4)*D_csf);
        %    S{ij} = f_wm*S_WM+f_gm*S_GM+f_csf*S_CSF;
        %    Ssim(index:index+ndirections(ij)-1) = S{ij};
        
        %
        index = index+ndirections(ij);
    end
    
else
    N = length(data.bvals);
    bmat{1} = zeros(N,4);
    bmat{1}(:,1:3) = data.bvecs;
    bmat{1}(:,4) = data.bvals;
    
    for i=1:length(phi)
        anglesFi = [phi(i), theta(i)]*(180/pi); % in degrees
        HRKernel(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K);
    end
    for l=1:length(isoDs)
        HRKernel(:,length(phi)+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end
    
    for i=1:length(phi_LR)
        anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
        LRKernel(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K);
    end
    for l=1:length(isoDs)
        LRKernel(:,length(phi_LR)+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end
end

end

% Private function assembling the NODDI based kernels.
function [bmat,LRKernel,HRKernel,super_scheme] = mDRLMT_MakeNODDIKernel_multicomp(data,nreconstruction_vertices,noddi_values,isoDs,shell_data)
if(isempty(which('SynthMeasWatsonSHStickTortIsoV_B0')))
    error('Cannot find the NODDI toolbox. Please, add it to the MATLAB path');
end

shells = single(unique(int16(round(data.bvals)/1)*1)); % automatic detection of number of shells.
% This shell splitting works only for b-value spaced more than 100. To be fixed for other datasets.
ndirections = zeros(length(shells),1);
for ij=1:length(shells)
    ndirections(ij) = sum(abs(data.bvals-shells(ij))<1);
end

bvals = data.bvals';
save('temp_bvals.bval','bvals','-ascii');
bvecs = data.bvecs';
% bvecs(:,1) = -bvecs(:,1);
% bvecs(:,3) = -bvecs(:,3);
save('temp_bvecs.bvec','bvecs','-ascii');
clear bvals bvecs
NODDI_protocol = FSL2Protocol('temp_bvals.bval','temp_bvecs.bvec');
delete('temp_bvals.bval');
delete('temp_bvecs.bvec');

lr_nreconstruction_vertices = min(length(data.bvals),90);

if(shell_data > 0)
    bmat = cell(length(shells),1); % bmat = [bvecs bvals]
    Kernel = cell(length(shells),1);
    Kernel_LR = cell(length(shells),1);
else
    bmat{1} = [data.bvecs data.bvals];
    Kernel = cell(1,1);
    Kernel_LR = cell(1,1);
end

% nreconstruction_vertices = 300; % 20/12/2017
nlr_vert = 0;
nhr_vert = sum(nreconstruction_vertices);

for ij=1:length(nreconstruction_vertices)
    super_scheme{ij} = gen_scheme(nreconstruction_vertices(ij),4); % the reconstruction scheme. Change 300 to any number
    [phi{ij}, theta{ij}] = cart2sph(super_scheme{ij}.vert(:,1),super_scheme{ij}.vert(:,2),super_scheme{ij}.vert(:,3)); % polar decomposition
    lr_scheme{ij} = gen_scheme(min(length(data.bvals),lr_nreconstruction_vertices(ij)),4);
    [phi_LR{ij}, theta_LR{ij}] = cart2sph(lr_scheme{ij}.vert(:,1),lr_scheme{ij}.vert(:,2),lr_scheme{ij}.vert(:,3));
    
    nlr_vert = nlr_vert + size(lr_scheme{ij}.vert,1);
end

HRKernel = zeros(sum(ndirections),nhr_vert+length(isoDs));
LRKernel = zeros(sum(ndirections),nlr_vert+length(isoDs));
% S = cell(length(shells),1);
% Ssim = zeros(sum(ndirections),1); % Just a simulated signal for internal testing

% noddi_values = [1 1.7E-9 3.5 0 3E-9 1];% x is the list of model parameters in SI units:
% % x(1) is the volume fraction of the intracellular space.
% % x(2) is the free diffusivity of the material inside and outside the cylinders.
% % x(3) is the concentration parameter of the Watson's distribution.
% % x(4) is the volume fraction of the isotropic compartment.
% % x(5) is the diffusivity of the isotropic compartment.
% % x(6) is the measurement at b=0.;

if(shell_data > 0)
    
    bvals = zeros(sum(ndirections),1);
    
    index = 1;
    
    for ij=1:length(shells)
        bmat{ij} = zeros(ndirections(ij),4);
        bmat{ij}(:,4) = shells(ij);
        bmat{ij}(:,1:3) = data.bvecs(index:index+ndirections(ij)-1,:);
        
        % Here the deconvolution dictionary is actually built
        
        % ANISOTROPIC PART
        hr_index_columns = 0;
        lr_index_columns = 0;
        
        Kernel{ij} = zeros(ndirections(ij),nhr_vert);
        Kernel_LR{ij} = zeros(ndirections(ij),nlr_vert);
        for aniso_comp = 1:length(nreconstruction_vertices)
            % HR
            for i=1:length(phi{aniso_comp})
                fibredir = super_scheme{aniso_comp}.vert(i,:)';
                E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
                Kernel{ij}(:,i+hr_index_columns) = E;
            end
            
            % LR
            for i=1:length(phi_LR{aniso_comp})
                fibredir = lr_scheme{aniso_comp}.vert(i,:)';
                E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
                Kernel_LR{ij}(:,i+lr_index_columns) = E;
            end
            
            hr_index_columns = hr_index_columns + size(super_scheme{aniso_comp}.vert,1);
            lr_index_columns = hr_index_columns + size(lr_scheme{aniso_comp}.vert,1);
        end
        
        % ISOTROPIC PART
        % HR
        for l=1:length(isoDs)
            Kernel{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end
        
        % LR
        for l=1:length(isoDs)
            Kernel_LR{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end
        
        % Build a linear dictionary with the subparts
        bvals(index:index+ndirections(ij)-1) = bmat{ij}(:,4);
        HRKernel(index:index+ndirections(ij)-1,:) = Kernel{ij};
        LRKernel(index:index+ndirections(ij)-1,:) = Kernel_LR{ij};
        
        index = index+ndirections(ij);
    end
    
else
    N = length(data.bvals);
    bmat{1} = zeros(N,4);
    bmat{1}(:,1:3) = data.bvecs;
    bmat{1}(:,4) = data.bvals;
    
    % ANISOTROPIC PART
    hr_index_columns = 0;
    lr_index_columns = 0;
    for aniso_comp = 1:length(nreconstruction_vertices)
        
        % HR
        for i=1:length(phi{aniso_comp})
            fibredir = super_scheme{aniso_comp}.vert(i,:)';
            E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
            HRKernel(:,i+hr_index_columns) = E;
        end
        %LR
        for i=1:length(phi_LR{aniso_comp})
            fibredir = lr_scheme{aniso_comp}.vert(i,:)';
            E = SynthMeasWatsonSHStickTortIsoV_B0(noddi_values{aniso_comp}, NODDI_protocol, fibredir);
            LRKernel(:,i+lr_index_columns) = E;
        end
        
        hr_index_columns = hr_index_columns + size(super_scheme{aniso_comp}.vert,1);
        lr_index_columns = lr_index_columns + size(lr_scheme{aniso_comp}.vert,1);
    end
    % ISOTROPIC PART
    % HR
    for l=1:length(isoDs)
        HRKernel(:,hr_index_columns+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end
    % LR
    for l=1:length(isoDs)
        LRKernel(:,lr_index_columns+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end
end

end

% The following functions have been included from "hardi_tools", which is
% available as follows:
% Original authors:
%% Project:   High Angular Resolution Diffusion Imaging Tools
% Available at: https://www.neuroimagen.es/webs/hardi_tools/
% Function to create a simulated signal from the multi-tensor diffusion model.
% Rician noise can be added.
%
%   Language:  MATLAB(R)
%   Author:  Erick Canales-Rodríguez, Lester Melie-García, Yasser Iturria-Medina, Yasser Alemán-Gómez
%   Date: 2013, Version: 1.2
%
% See also test_DSI_example, test_DOT_example, test_QBI_example,
% test_DOT_R1_example, test_DOT_R2_vs_CSA_QBI_example.

% Modifications have been performed where needed for compatibility (from
% A. De Luca)
function [S, D] = create_signal_multi_tensor (ang, f, Eigenvalues, b, grad, S0, SNR, add_noise)
% -Normalizing the gradient vector and then transforming the b-value.
% -This part is only for non-unitary gradient vectors
% Transform_b_value_and_grad = repmat(sqrt(diag(grad*grad')+eps), [1 3]);
% grad = grad./Transform_b_value_and_grad;
% b = b.*(Transform_b_value_and_grad(:,1)).^2;

A = diag(Eigenvalues);

S = 0;
Nfibers = length(f);
f = f/sum(f);
for i = 1:Nfibers
    phi(i) = ang(i, 1);
    theta(i) = ang(i, 2);
    R = RotMatrix(phi(i),theta(i));
    D = R*A*R';
    S = S + f(i)*exp(-b.*diag(grad*D*grad'));
end
S = S0*S;

sigma = S0/SNR;

standar_deviation = sigma.*(ones(length(grad),1));
med = zeros(length(grad),1);

er1 = normrnd(med, standar_deviation);
er2 = normrnd(med, standar_deviation);

if add_noise == 1
    S = sqrt((S + er1).^2 + er2.^2);
end

return
end

% Extension to DKI - A. De Luca
function [S, D] = create_signal_multi_tensor_dki (ang, f, Eigenvalues, b, grad, S0, SNR, add_noise, K)
% -Normalizing the gradient vector and then transforming the b-value.
% -This part is only for non-unitary gradient vectors
% Transform_b_value_and_grad = repmat(sqrt(diag(grad*grad')+eps), [1 3]);
% grad = grad./Transform_b_value_and_grad;
% b = b.*(Transform_b_value_and_grad(:,1)).^2;

A = diag(Eigenvalues);

S = 0;
Nfibers = length(f);
f = f/sum(f);
b2D2K = 1/6*b.^2.*mean(Eigenvalues.^2)*K;
for i = 1:Nfibers
    phi(i) = ang(i, 1);
    theta(i) = ang(i, 2);
    R = RotMatrix(phi(i),theta(i));
    D = R*A*R';
    S = S + f(i)*exp(-b.*diag(grad*D*grad')+b2D2K);
end
S = S0*S;

sigma = S0/SNR;

standar_deviation = sigma.*(ones(length(grad),1));
med = zeros(length(grad),1);

er1 = normrnd(med, standar_deviation);
er2 = normrnd(med, standar_deviation);

if add_noise == 1
    S = sqrt((S + er1).^2 + er2.^2);
end
end

% --- private funtions -----------
function R = RotMatrix(phi,theta)

c = pi/180;
phi = phi*c;
theta = theta*c;

Rz = [ cos(phi)  -sin(phi)  0
    sin(phi)   cos(phi)  0
    0         0      1];


Ry = [cos(theta)   0   sin(theta)
    0        1         0
    -sin(theta)   0   cos(theta)];

R =  Rz*Ry;
return
end

function scheme = gen_scheme(N, lmax)

% function scheme = gen_scheme(N, lmax)
%
% Generate a set of orientations in the required format, along
% with the corresponding SH transform information up to
% harmonic order 'lmax'.
%
% If N is a string, it will attempt to load the specified
% file.
%
% If N is a number, a scheme with the specified number of
% directions will be generated using the equidistribute.m
% script (note that these are not perfectly uniformly
% distributed).
%
% If N is a nx3 matrix, it will assume that each row provides
% an [ x y z ] vector pointing along the desired direction.
%

if ischar(N)
    N = load(N);
end

if size(N,1) == 1 & size(N,2) == 1
    P = c2s(equidistribute(N));
elseif size(N,2) >= 3
    n = sqrt(sum(N(:,1:3).^2,2));
    k = find(n);
    X = N(k,1:3)./repmat(n(k),1,3);
    P = c2s(X);
else
    P = N;
end

scheme.el = P(:,1);
scheme.az = P(:,2);

scheme.sh = [];
scheme.lmax = lmax;

for l = 0:2:lmax
    scheme.sh = [ scheme.sh eval_SH(l, scheme.el, scheme.az)' ];
end

scheme.vert= s2c([ scheme.el scheme.az 1+0*scheme.az]);
scheme.mesh = convhulln(scheme.vert);
end

function s = eval_ALP(l, el)

% syntax: s = eval_ALP(l, el)
%
% evaluates the Associated Legendre Polynomial at elevations 'el'
% for harmonic order l.

s = legendre(l, cos(el'));
for m = 0:l
    s(m+1,:) = s(m+1,:).*sqrt((2*l+1)*factorial(l-m) / ((4*pi)*factorial(l+m)));
end

if l
    s = [ s(end:-1:2,:); s ];
end
end

function s = eval_SH(l, el, az)

% syntax: s = eval_SH(l, el, az)
%
% Evaluates the lth order spherical harmonic coefficients at
% positions [ el az ].

s = ones(size(az,1),1);

if l > 0
    s = [ sqrt(2)*sin(az*(l:-1:1)) s sqrt(2)*cos(az*(1:l)) ];
end

s = eval_ALP(l, el).*s';
end

function X = equidistribute(N)
% X = equidistribute(N)
% uses the formula in [saff:dmp:1997] to generate equidistributed
% points on the sphere.
% INPUT: N is the number of points, default=12
%
% OUTPUT: X is the set of points as a Nxd matrix, one point per row
%         if no output is specified, the points are plotted using scatter3.
% REFERENCE
% @Article{saff:dmp:1997,
%  author =          {Saff, E. B. and Kuijlaars, A. B. J.},
%  title =          {Distributing Many Points on a Sphere},
%  journal =          {Math. Intell.},
%  year =          {1997},
%  volume =          {19},
%  number =          {1},
%  pages =          {5--11},
%}


X = zeros(N,3);

for k=1:N
    h = -1 + 2*(k-1)/(N-1);
    theta(k) = acos(h);
    if k==1 | k==N
        phi(k) = 0;
    else
        phi(k) = mod(phi(k-1) + 3.6/sqrt(N*(1-h^2)),2*pi);
    end;
    X(k,:) = [ cos(phi(k))*sin(theta(k)), ...
        sin(phi(k))*sin(theta(k)), ...
        cos(theta(k)) ];
end;

%if nargout == 0
%Z = zeros(size(X));
%[SX,SY,SZ] = sphere;
%scatter3(X(:,1),X(:,2),X(:,3),'filled')
%hold on
%sph2 = surf(SX,SY,SZ);
%set(sph2, 'FaceColor', [ 1 1 0 ]);
%axis vis3d
%end
end

% Taken from "hardi tools" of Erick Canales-Rodriguez
function fODF = ADT_deconv_RLdamp_1D_noEP(Signal, Kernel, Niter,nu)
% fODF: purely anisotropic part

fODF0 = ones(size(Kernel,2),1);

fODF = fODF0/sum(fODF0);
KernelT = Kernel'; % only one time

fzero = 1e-06;
KernelS = KernelT*Signal;
V = 8;

mu = max(0, 1 - 4*std(Signal) );
nuV = nu^V;
last_odf = fODF;
my_eps = 1e-4*max(Signal);
for i = 1:Niter
    % --- approach: Flavio Dell’Acqua
    Dem = KernelT*(Kernel*fODF);
    fODFV = fODF.^V;
    Rk = 1 - (fODFV)./(fODFV + nuV);
    Uk = 1 - mu*Rk;
    RL_factor = 1 + Uk.*( (KernelS-Dem)./(Dem + eps) );
    fODF = max(fzero, fODF.*RL_factor); % positivity
    if(sum(abs(fODF-last_odf)) <= my_eps)
        break
    end
    last_odf = fODF;
end


end

% My implementation of the classic Richardson Lucy, inspired by the above.
% Thanks to Flavio Dell'Acqua for help in drafting these methods.
function fODF = RichardsonLucy(Signal, Kernel, Niter)

fODF0 = ones(size(Kernel,2),1);

fODF = fODF0/sum(fODF0);
KernelT = Kernel'; % only one time

KernelS = KernelT*Signal;

last_odf = fODF;
my_eps = 1e-4*max(Signal);
for i = 1:Niter
    Dem = KernelT*Kernel*fODF;
    fODF = fODF.*(KernelS./(Dem+eps));
    if(sum(abs(fODF-last_odf)) <= my_eps)
        break
    end
    last_odf = fODF;
end


end

function nn_val = get_drl_nn_heuristic(bh,bvalue,tgt_diffusivity)
if(nargin < 3)
    tgt_diffusivity = 0.7e-3;
end
S = exp(-bvalue*tgt_diffusivity)*ones(size(bh,1),1);
K = RichardsonLucy(S,bh,200);
nn_val = max(K)*2;
end


% Taken from the MRtrix package
function tracks = read_mrtrix_tracks (filename)

% function: tracks = read_mrtrix_tracks (filename)
%
% returns a structure containing the header information and data for the MRtrix 
% format track file 'filename' (i.e. files with the extension '.tck'). 
% The track data will be stored as a cell array in the 'data' field of the
% return variable.

f = fopen (filename, 'r');
assert(f ~= -1, 'error opening %s', filename);
L = fgetl(f);
if ~strncmp(L, 'mrtrix tracks', 13)
  fclose(f);
  error('%s is not in MRtrix format', filename);
end

tracks = struct();

while 1
  L = fgetl(f);
  if ~ischar(L), break, end;
  L = strtrim(L);
  if strcmp(L, 'END'), break, end;
  d = strfind (L,':');
  if isempty(d)
    disp (['invalid line in header: ''' L ''' - ignored']);
  else
    key = lower(strtrim(L(1:d(1)-1)));
    value = strtrim(L(d(1)+1:end));
    if strcmp(key, 'file')
      file = value;
    elseif strcmp(key, 'datatype')
      tracks.datatype = value;
    else 
      tracks = add_field (tracks, key, value);
    end
  end
end
fclose(f);

assert(exist('file') && isfield(tracks, 'datatype'), ...
  'critical entries missing in header - aborting');

[ file, offset ] = strtok(file);
assert(strcmp(file, '.'), ...
  'unexpected file entry (should be set to current ''.'') - aborting');

assert(~isempty(offset), 'no offset specified - aborting');
offset = str2num(char(offset));

datatype = lower(tracks.datatype);
byteorder = datatype(end-1:end);

if strcmp(byteorder, 'le')
  f = fopen (filename, 'r', 'l');
  datatype = datatype(1:end-2);
elseif strcmp(byteorder, 'be')
  f = fopen (filename, 'r', 'b');
  datatype = datatype(1:end-2);
else
  error('unexpected data type - aborting');
end

assert(f ~= -1, 'error opening %s', filename);

fseek (f, offset, -1);
data = fread(f, inf, datatype);
fclose (f);

N = floor(prod(size(data))/3);
data = reshape (data, 3, N)';
k = find (~isfinite(data(:,1)));

tracks.data = {};
pk = 1;
for n = 1:(prod(size(k))-1)
  tracks.data{end+1} = data(pk:(k(n)-1),:);
  pk = k(n)+1;
end
end

% Taken from MRtrix
function image = add_field (image, key, value)
  if isfield (image, key)
    previous = getfield (image, key);
    if iscell (previous)
      image = setfield (image, key, [ previous {value} ]);
    else
      image = setfield (image, key, { previous, value });
    end
  else
    image = setfield (image, key, value);
  end
end

% taken from MRtrix
function write_mrtrix_tracks (tracks, filename)

% function: write_mrtrix_tracks (tracks, filename)
%
% writes the track data stored as a cell array in the 'data' field of the
% tracks variable to the MRtrix format track file 'filename'. All other fields
% of the tracks variable will be written as text entries in the header, and are
% expected to supplied as character arrays.

assert(isfield(tracks, 'data'), ...
  'input tracks variable does not contain required ''data'' field');

assert(iscell(tracks.data), ...
  'input tracks.data variable should be a cell array');

f = fopen (filename, 'w', 'ieee-le');
assert(f ~= -1, 'error opening %s', filename);

fprintf (f, 'mrtrix tracks\ndatatype: Float32LE\ncount: %d\n', prod(size(tracks.data)));
names = fieldnames(tracks);
for i=1:size(names)
  if strcmpi (names{i}, 'data'), continue; end
  if strcmpi (names{i}, 'count'), continue; end
  if strcmpi (names{i}, 'datatype'), continue; end
  fprintf (f, '%s: %s\n', names{i}, getfield(tracks, names{i}));
end
data_offset = ftell (f) + 20;
fprintf (f, 'file: . %d\nEND\n', data_offset);

fwrite (f, zeros(data_offset-ftell(f),1), 'uint8');
for i = 1:prod(size(tracks.data))
  fwrite (f, tracks.data{i}', 'float32');
  fwrite (f, [ nan nan nan ], 'float32');
end

fwrite (f, [ inf inf inf ], 'float32');
fclose (f);
end

% Function to convert from MAT 2 VTK. Original implementation in
% ExploreDTI, adapted and fixed from A. De Luca
function [mx,my,mz] = E_DTI_Convert_tracts_mat_2_vtk_lines(filename, filename_out)

try
    warning off all
    load(filename,'Tracts','TractMask')
    warning on all
catch
    return;
end

if ~exist('Tracts','var') || ~iscell(Tracts)
    return;
end

Tract_dir = E_DTI_Get_tract_dir(Tracts);

for i=1:length(Tract_dir)
    Tract_dir{i} = abs(Tract_dir{i});    
end

num_tracts = length(Tracts);
TL = zeros(num_tracts,1);

for i=1:num_tracts
    TL(i)=size(Tracts{i},1);
end

CTL = [0; cumsum(TL)];
CTL_c = [0; cumsum(TL+1)];

Total_P = sum(TL);
Points = repmat(single(0),[Total_P 3]);
% Colors = repmat(single(1),[Total_P 4]);
Colors = repmat(single(1),[Total_P 3]);

for i=1:num_tracts
    Points(CTL(i)+1:CTL(i+1),:) = single(Tracts{i});
    Points(CTL(i)+1:CTL(i+1),[2 1 3]) = single(Tracts{i}); % mod ADL
    Points(CTL(i)+1:CTL(i+1),1) = size(TractMask,2) - Points(CTL(i)+1:CTL(i+1),1); % mod ADL
    Points(CTL(i)+1:CTL(i+1),2) = size(TractMask,1) - Points(CTL(i)+1:CTL(i+1),2); % mod ADL

    Colors(CTL(i)+1:CTL(i+1),[2 1 3]) = single(Tract_dir{i});
end

mx = mean(Points(:,1));
my = mean(Points(:,2))+40;
mz = mean(Points(:,3))-10;

Points(:,1) = Points(:,1) - mx;
Points(:,2) = Points(:,2) - my;
Points(:,3) = Points(:,3) - mz;

conn = repmat(single(0),[1 num_tracts + Total_P]);

for i=1:num_tracts
    conn(CTL_c(i)+1) = TL(i);
    conn(CTL_c(i)+2:CTL_c(i+1)) = (CTL(i):CTL(i+1)-1);
end

if nargin==1
    [PATHSTR,NAME,EXT] = fileparts(filename);
    save_name = [PATHSTR filesep NAME '.vtk'];
else
    save_name = filename_out;
end

try
    
    fid = fopen(save_name,'w+t');
    
    fprintf(fid, '%s\n','# vtk DataFile Version 3.0');
    fprintf(fid, '%s\n','vtk output');
    fprintf(fid, '%s\n','ASCII');
    fprintf(fid, '%s\n','DATASET POLYDATA');
    fprintf(fid, ['POINTS ' '%d' ' float' '\n'], Total_P);
    fprintf(fid,'%f %f %f\n',Points');
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid, '\n');
    fprintf(fid, ['LINES ' '%d %d\n'], [num_tracts num_tracts + Total_P]);
    fprintf(fid,'%d ',conn);
    fprintf(fid, '\n\n');
    fprintf(fid, ['POINT_DATA ' '%d' '\n'], Total_P);
    % fprintf(fid, '%s\n','COLOR_SCALARS Vertex_colors 4');
    % fprintf(fid,'%f %f %f %f\n',Colors');
    fprintf(fid, '%s\n','COLOR_SCALARS rgb_colors 3');
    fprintf(fid,'%f %f %f\n',Colors');
    fclose(fid);
    
catch
    disp('Problems writing data... check permissions...')
    return;
end

end

% Helper function from EDTI
function Tract_dir = E_DTI_Get_tract_dir(Tracts)

if length(Tracts)==0
    Tract_dir = cell(1,0);
end

for i=1:length(Tracts)
    
    Tract_dir{i} = Tracts{i};
    
    if size(Tracts{i},1)>2
        dummy = Tract_dir{i}(2:end,:)-Tract_dir{i}(1:end-1,:);
        Tract_dir{i}(2:end-1,:) = (dummy(1:end-1,:)+dummy(2:end,:))/2;
        Tract_dir{i}(1,:) = dummy(1,:);
        Tract_dir{i}(end,:) = dummy(end,:);
    else
        dummy = Tract_dir{i}(2,:)-Tract_dir{i}(1,:);
        Tract_dir{i}(1,:) = dummy;
        Tract_dir{i}(2,:) = dummy;
    end
    
    
    NT = sqrt(sum(Tract_dir{i}.^2,2));
    
    Tract_dir{i} = Tract_dir{i}./[NT NT NT];
    
    
end
end

% A simple VTK read function
function [points_list,lines] = VTKImport_test(file_in)

fid = fopen(file_in,'rb');

str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
nvert = sscanf(str,'%*s %d %*s', 1);
sp = ftell(fid);
% read vertices
A = single(fread(fid,nvert*3,'float','b'));
mp = ftell(fid);

str = fgets(fid);
str = fgets(fid);
nlines = sscanf(str,'%*s %d %d');
lines = cell(1,1);
for ij=1:nlines
    K = fread(fid,1,'int','b');
    lines(ij) = {fread(fid,K,'*int','b')};
end
% B = fread(fid,[nlines(1) 4],'*uint8','b');
pstr = fgets(fid);
mmp = ftell(fid);

fclose(fid);

A = reshape(A,3,nvert);
points_list = A';
end
