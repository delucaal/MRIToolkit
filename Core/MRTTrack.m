%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

% MRTTrack (former Spherical deconvolution) class - including GRL and mFOD
% A. De Luca - UMC Utrecht - 28/10/2019 - alberto@isi.uu.nl -
% alberto.deluca.06@gmail.com
% First version: 28/10/2019

classdef MRTTrack < handle
    properties
        data;
        n_isotropic;
        n_anisotropic;
        deconv_method;
        inner_shells_weight;
        shell_weight;
        rf_indices_lr;
        rf_indices_hr;
        rf_models;
        L2LSQ_reg = 0.1;
        LRKernel;
        HRKernel;
        weighted_LRKernel;
        weighted_HRKernel;
        nDirections;
        NN_L, NN_H;
        data_size;
        % for robust deconvolution (outlier rejection)
        max_outlier_rejection = 0.1; % 10%
        robust_kappa = 3;
        robust_deconv  = 0;% whether to use robust_deconv
        Lmax_FOD = 8;
        drl_iters = 200;
        stdReconDirections;
    end

    methods

        function obj = MRTTrack(varargin)
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
            gs = gen_scheme(300,4);
            obj.stdReconDirections = gs.vert;

            if(length(varargin) > 0 && iscell(varargin{1}) > 0)
                varargin = varargin{1};
            end
            InputPairs = ParseInputKeys(varargin);
            for ij=1:size(InputPairs,1)
                if(strcmpi(InputPairs{ij,1},'data'))
                    obj.data = InputPairs{ij,2};
                    obj.data_size = size(obj.data.img);
                end
            end
        end

        function obj = setInnerShellWeighting(obj,weight)
            %       % How the inner shells are weighted in the deconvolution (0-1)
            obj.inner_shells_weight = weight;
            obj.shell_weight = zeros(size(obj.data.bvals));
            all_shells = (unique(obj.data.bvals));
            for ij=1:length(all_shells)
                this_shell = abs(obj.data.bvals-all_shells(ij)) < 100;
                obj.shell_weight(this_shell) = ij;
            end

            obj.shell_weight(obj.shell_weight < length(all_shells)) = obj.inner_shells_weight;
            obj.shell_weight(obj.shell_weight == length(all_shells)) = 1;
            obj.weighted_LRKernel = obj.LRKernel;
            obj.weighted_HRKernel = obj.HRKernel;

            % Weighted kernels for the lower shells
            for ij=1:size(obj.LRKernel,2)
                obj.weighted_LRKernel(:,ij) = obj.LRKernel(:,ij).*obj.shell_weight;
            end

            for ij=1:size(obj.HRKernel,2)
                obj.weighted_HRKernel(:,ij) = obj.HRKernel(:,ij).*obj.shell_weight;
            end

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
            [~,lLRKernel,lHRKernel] = mDRLMT_MakeDKIKernel_multicomp(obj.data,10,[1.7e-3 0.5e-3 0.5e-3],0,0,D,1);
            if(~isempty(obj.LRKernel))
                obj.LRKernel(:,IndexLR) = single(lLRKernel(:,end));
                obj.HRKernel(:,IndexHR) = single(lHRKernel(:,end));
            else
                obj.LRKernel = single(lLRKernel(:,end));
                obj.HRKernel = single(lHRKernel(:,end));
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

        function obj = AddAnisotropicRF_DKI(obj,EigenValues,MeanKurtosis,Offset)
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

            if(exist('Offset','var') < 1)
                Offset = 0;
            end
            
            pass_data = obj.data;
            pass_data.bvecs(:,3) = -pass_data.bvecs(:,3);
            obj.n_anisotropic = obj.n_anisotropic+1;
            [~,lLRKernel,lHRKernel] = mDRLMT_MakeDKIKernel_multicomp(pass_data,obj.nDirections,EigenValues,MeanKurtosis,Offset,[],0);
            IndexHR = size(obj.HRKernel,2)+1:size(obj.HRKernel,2)+obj.nDirections;
            IndexLR = size(obj.LRKernel,2)+1:size(obj.LRKernel,2)+size(lLRKernel,2);
            if(~isempty(obj.LRKernel))
                obj.LRKernel(:,IndexLR) = single(lLRKernel);
                obj.HRKernel(:,IndexHR) = single(lHRKernel);
            else
                obj.LRKernel = single(lLRKernel);
                obj.HRKernel = single(lHRKernel);
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
                obj.LRKernel(:,IndexLR) = single(lLRKernel);
                obj.HRKernel(:,IndexHR) = single(lHRKernel);
            else
                obj.LRKernel = single(lLRKernel);
                obj.HRKernel = single(lHRKernel);
            end
            obj.rf_indices_lr(end+1) = {IndexLR};
            obj.rf_indices_hr(end+1) = {IndexHR};
            obj.rf_models(end+1) = {'NODDI'};
        end

        function obj = AutomaticDRLDamping(obj)
            % Compute the automatic damping for the dRL method. Only needed if
            % the deconvolution method is dRL.
            TGT_ADC = 0.7e-3;
            %             if(~isempty(obj.rf_models))
            %                for rf_id=1:length(obj.rf_models)
            %                   if(strcmp(obj.rf_models{rf_id},'ADC') > 0)
            %                       TGT_ADC =
            %                   end
            %                end
            %             end
            obj.NN_L = get_drl_nn_heuristic(obj.LRKernel,max(obj.data.bvals),TGT_ADC);
            obj.NN_H = get_drl_nn_heuristic(obj.HRKernel,max(obj.data.bvals),TGT_ADC);
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

            obj.setInnerShellWeighting(obj.inner_shells_weight);
            %             weighted_LRKernel = obj.LRKernel;
            %             weighted_HRKernel = obj.HRKernel;

            fractions = zeros([sm NC+ANC],'single');
            if(low_mem_mode == 0)
                WM_fod = zeros([sm nreconstruction_vertices],'single');
            else
                WM_fod = sparse([sm nreconstruction_vertices],'single');
            end
            RSS = zeros([sm 1],'single');
%             S0 = zeros([sm 1],'single');

            N = sm;
            op_e2 = optimset('TolX',1e-2);
            op = optimset('TolX',1e-6);
            max_discarded_vols = round(obj.max_outlier_rejection*st);

            [~,DeconvMethodCode] = SphericalDeconvolution.isSupportedMethod(obj.deconv_method); % -1 = failure; 1 = LSQNONNEG; 2 = DW_RegularizedDeconv; 3 = dRL
            if(DeconvMethodCode == -1)
                warning('Unsupported deconvolution method.');
                return;
            end

            TheTol = 1e-3;
            tic

            min_bval = min(obj.data.bvals);

            if(NC > 0 || ANC > 1) % there are multiple isotropic or anisotropic components (otherwise run dRL)
                fODF0_LR = zeros(nreconstruction_vertices_lr,1,'single');
                fODF0_LR = fODF0_LR/sum(fODF0_LR);
                fODF0_HR = zeros(nreconstruction_vertices,1,'single');
                fODF0_HR = fODF0_HR/sum(fODF0_HR);

                LLRKernel = obj.weighted_LRKernel(:,1:end-NC);
                ULRKernel = obj.LRKernel(:,1:end-NC);
                LLRKernelT = LLRKernel';
                LLKTK = LLRKernelT*LLRKernel;
                LHRKernel = obj.weighted_HRKernel(:,1:end-NC);
                UHRKernel = obj.HRKernel(:,1:end-NC);
                LHRKernelT = LHRKernel';
                LHKTK = LHRKernelT*LHRKernel;
                LLRKernelIso = obj.weighted_LRKernel(:,end-NC+1:end);
                LHRKernelIso = obj.weighted_LRKernel(:,end-NC+1:end);
                ULRKernelIso = obj.LRKernel(:,end-NC+1:end);
                UHRKernelIso = obj.LRKernel(:,end-NC+1:end);

                max_voxels_per_worker = 3e4;
                mask_vector = find(obj.data.mask > 0);
                nruns = ceil(length(mask_vector)/max_voxels_per_worker);
                for run=1:nruns
                    disp(['Run ' num2str(run) ' of ' num2str(nruns)])
                    indices = (1+(run-1)*max_voxels_per_worker):min(run*max_voxels_per_worker,...
                        length(mask_vector));
                    Nindices = length(indices);
                    obj_copy = MRTTrack('data',obj.data);
                    props = properties(obj);
                    for k = 1:length(props)
                        if isprop(obj_copy,props{k})
                           obj_copy.(props{k}) = obj.(props{k});
                        end
                    end
                    obj_copy.data.img = obj_copy.data.img(:,mask_vector(indices));
                    obj_copy.data.mask = obj_copy.data.mask(mask_vector(indices));
                    l_fractions = zeros([Nindices size(fractions,2)],'single');
                    l_WM_fod = zeros([Nindices size(WM_fod,2)],'single');
                    l_RSS = zeros([Nindices 1],'single');

                    parfor (x=1:Nindices,nof_workers)
%                     for x=1:N
                        YLR = [ULRKernel(:,1) ULRKernelIso];
                        YHR = [UHRKernel(:,1) UHRKernelIso];
                        %                     fODFC = [];
                        fODF = [];
%                         if(obj_copy.data.mask(x) < 1)
%                             continue
%                         end

                        Stot = obj_copy.data.img(:,x);
                        fit_data = true(size(Stot));

                        % This normalization assumes there is some b=0s/mm2 data. It is not
                        % essential to do it as long as fractions is normalized at the end
                        NormFactor = mean(Stot(obj_copy.data.bvals<=min_bval));
                        %                         S0(x) = NormFactor;
                        Stot = Stot/NormFactor;

                        Stot_w = Stot.*obj_copy.shell_weight; % Weight the lower shells

                        piso = zeros(NC,1);
                        p_old = Inf; % store the fractions at previous iter

                        % The following loop will iterate WM-FOD estimation (with mDRL) and
                        % fractions guess (with LSQNONNEG)
                        fODFC = fODF0_LR;
                        for iter=1:50 % 50 = max number of iterations but usually exits way before

                            DS = max(Stot_w(fit_data) - LLRKernelIso(fit_data,:)*piso,0); % Subtract GM and CSF contributions

                            if(DeconvMethodCode == 4)
                                %                         fODFC = mat_dRL(DS, weighted_LRKernel(:,1:end-NC), 200, obj.NN_L, 8);
                                fODFC = ADT_deconv_RLdamp_1D_noEP_init(fODF0_LR,DS, LLRKernel(fit_data,:), LLRKernelT(:,fit_data), LLKTK, 200, obj_copy.NN_L); %obj_copy.drl_iters
%                                 fODFC = ADT_deconv_RLdamp_1D_noEP(DS, LLRKernel(fit_data,:), obj_copy.drl_iters, obj_copy.NN_L); %obj_copy.drl_iters
                            elseif(DeconvMethodCode == 3)
                                %                         fODFC = mat_RL(DS, weighted_LRKernel(:,1:end-NC), 200);
                                fODFC = RichardsonLucy(DS, LLRKernel(fit_data,:), 200);
                                %                         elseif(DeconvMethodCode == 2)
                                %                             fODFC = DW_RegularizedDeconv(weighted_LRKernel(fit_data,1:end-NC),DS,op_e2,obj.L2LSQ_reg);
                            elseif(DeconvMethodCode == 2 || DeconvMethodCode == 1)
                                fODFC = lsqnonneg(double(LLRKernel(fit_data,:)),double(DS),op);
                                %fODFC = DW_RegularizedDeconv(obj.weighted_LRKernel(fit_data,1:end-NC),DS,op, obj.L2LSQ_reg);
                            elseif(DeconvMethodCode == 6)
                                fODFC = sph_deconv_tv_motor_all(DS, LLRKernel(fit_data,:), zeros(size(obj_copy.weighted_LRKernel,2)-NC,1), obj_copy.drl_iters, 1, 'SMF-SENSE-based');
                            elseif(DeconvMethodCode == 7)
                                fODFC = lasso(LLRKernel(fit_data,:),DS,'Alpha',0.5,'Lambda',obj_copy.L2LSQ_reg);
                            elseif(DeconvMethodCode == 8)
                                fODFC = MRT_Library.ConstrainedDeconvolution(LLRKernel(fit_data,:),DS,obj_copy.L2LSQ_reg);
                            end
                            % This line is quite tricky. It comes out of some trial and error
                            % but actually has the intent of eliminating 1) small contributions
                            % 2) flat - spread fods enforcing sparsity.
                            % Without this line the code DOESN'T work. (WM is over-estimated).
                            fODFC(fODFC < median(fODFC)) = 0;
                            %                         fODFC(fODFC < 0.01 | fODFC < 0.1*max(fODFC)) = 0;
                            % Build a dictionary to fit the complete signal (Stot)
                            YLR(:,1) = ULRKernel*fODFC; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                            M = max(YLR(:,1));
                            if(M > 0) % i.e. if the FOD is non-zero
                                YLR(:,1) = YLR(:,1)/M; % Normalize WM signal
                                %                               fODFC = fODFC / sum(fODFC);
                            end

                            p = lsqnonneg(double(YLR(fit_data,:)),double(Stot(fit_data))); % Compute the signal fractions
                            piso = p(end-NC+1:end);

                            % if nothing changed compared to previous iter, exit. (Tol may need to be
                            % adjusted)
                            if(sum(abs(p-p_old) < TheTol) == 3)
                                break
                            end
                            p_old = p;

                            if(obj_copy.robust_deconv == 1 && sum(fit_data == 0) < max_discarded_vols)
                                residuals = Stot_w - Y*p;
                                fit_data = abs(residuals - median(residuals)) < obj_copy.robust_kappa*1.4826*mad(residuals,1);
                            end

                        end

                        % NEW FINAL STEP 05/02/2018
                        DS = max(Stot_w(fit_data) - LHRKernelIso(fit_data,:)*piso,0); % Subtract GM and CSF contributions

                        if(DeconvMethodCode == 4)
                            %                     fODF = mat_dRL(DS, weighted_HRKernel(:,1:end-NC),200, obj.NN_H, 8);
                            fODF = ADT_deconv_RLdamp_1D_noEP_init(fODF0_HR,DS, LHRKernel(fit_data,:), LHRKernelT(:,fit_data), LHKTK, obj_copy.drl_iters, obj_copy.NN_H);
%                             fODF = ADT_deconv_RLdamp_1D_noEP(DS, LHRKernel(fit_data,:), obj_copy.drl_iters, obj_copy.NN_H);
                        elseif(DeconvMethodCode == 3)
                            %                     fODF = mat_RL(DS, weighted_HRKernel(:,1:end-NC), 200);
                            fODF = RichardsonLucy(DS, LHRKernel(fit_data,:), obj_copy.drl_iters);
                        elseif(DeconvMethodCode == 2)
                            fODF = DW_RegularizedDeconv(double(LHRKernel(fit_data,:)),double(DS),op, obj_copy.L2LSQ_reg);
                        elseif(DeconvMethodCode == 1)
                            fODF = lsqnonneg(double(LHRKernel(fit_data,:)),double(DS));
                        elseif(DeconvMethodCode == 6)
                            fODF = sph_deconv_tv_motor_all(DS, LHRKernel(fit_data,:), zeros(size(obj_copy.weighted_HRKernel,2)-NC,1), obj_copy.drl_iters, 1, 'SMF-SENSE-based');
                        elseif(DeconvMethodCode == 7)
                            fODF = lasso(LHRKernel(fit_data,:),DS,'Alpha',0.5,'Lambda',obj_copy.L2LSQ_reg);
                        elseif(DeconvMethodCode == 8)
                            fODF = MRT_Library.ConstrainedDeconvolution(LHRKernel(fit_data,:),DS,obj_copy.L2LSQ_reg);
                        end
                        fODFC = fODF;
                        fODFC(fODFC < median(fODFC)) = 0;
                        %                     fODFC(fODFC < 0.01 | fODFC < 0.1*max(fODFC)) = 0;
                        %                     fODF = fODFC;

                        YHR(:,1) = UHRKernel*fODFC; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                        if(obj_copy.n_anisotropic == 1)
                            M = max(YHR(:,1));
                            if(M > 0) % i.e. if the FOD is non-zero
                                %                             fODFC = fODFC / sum(fODFC);
                                YHR(:,1) = YHR(:,1)/M; % Normalize WM signal
                            end
                        else
                            YHR = [];
                            for kANC = 1:ANC
                                NewCol = UHRKernel(:,obj_copy.rf_indices_hr{kANC})*fODFC(obj_copy.rf_indices_hr{kANC});
                                M = max(NewCol);
                                if(M > 0)
                                    NewCol = NewCol/M;
                                end
                                YHR = [YHR NewCol];
                            end
                            YHR = [YHR UHRKernelIso]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                        end
                        p = lsqnonneg(double(YHR(fit_data,:)),double(Stot(fit_data)),op_e2); % Compute the signal fractions
                        l_RSS(x) = sum((Stot_w-YHR*p).^2);

                        l_fractions(x,:) = p;
                        l_WM_fod(x,:) = single(fODF);
                    end
                    WM_fod(mask_vector(indices),:) = l_WM_fod;
                    fractions(mask_vector(indices),:) = l_fractions;
                    RSS(mask_vector(indices)) = l_RSS;
                end
            else
                max_voxels_per_worker = 3e4;
                nruns = ceil(size(obj.data.img,2)/max_voxels_per_worker);
                for run=1:nruns
                    disp(['Run ' num2str(run) ' of ' num2str(nruns)])
                    %                 parfor (x=1:N,nof_workers)
                    indices = (1+(run-1)*max_voxels_per_worker):min(run*max_voxels_per_worker,...
                        size(obj.data.img,2));
                    Nindices = length(indices);
                    obj_copy = MRTTrack('data',obj.data);
                    props = properties(obj);
                    for k = 1:length(props)
                        if isprop(obj_copy,props{k})
                           obj_copy.(props{k}) = obj.(props{k});
                        end
                    end
                    obj_copy.data.img = obj_copy.data.img(:,indices);
                    obj_copy.data.mask = obj_copy.data.mask(indices);
                    l_fractions = zeros([Nindices size(fractions,2)],'single');
                    l_WM_fod = zeros([Nindices size(WM_fod,2)],'single');
                    l_RSS = zeros([Nindices 1],'single');

                    parfor (x=1:Nindices,nof_workers)
                        %     if(mod(x,progressStepSize) == 0)
                        %         ppm.increment();
                        %     end
                        if(obj_copy.data.mask(x) < 1)
                            continue
                        end

                        Stot = double(obj_copy.data.img(:,x));

                        % This normalization assumes there is some b=0s/mm2 data. It is not
                        % essential to do it as long as fractions is normalized at the end
                        NormFactor = mean(Stot(obj_copy.data.bvals<=min_bval));
                        %                     S0(x) = NormFactor;
                        Stot = Stot/NormFactor;

                        l_fractions(x,1) = 1;
                        FOD = ADT_deconv_RLdamp_1D_noEP(Stot, obj_copy.weighted_HRKernel,obj_copy.drl_iters, obj_copy.NN_H/100);
                        l_WM_fod(x,:) = single(FOD);
                        l_RSS(x) = sum((Stot-obj_copy.weighted_HRKernel*FOD).^2);
                    end
                    WM_fod(indices,:) = l_WM_fod;
                    fractions(indices,:) = l_fractions;
                    RSS(indices) = l_RSS;

                end
            end
            % Reshape data
            toc

            fsum = sum(fractions,2);
            for ij=1:size(fractions,2)
                fractions(:,ij) = fractions(:,ij) ./ (fsum+eps);
            end

            output.fractions = reshape(fractions,[siz(1:3),size(fractions,2)]);

            %             WM_fod_max = max(WM_fod,[],2);
            %             WM_fod_normalized = WM_fod;
            %             WM_fod_val = mean(WM_fod_max(fractions(:,1) > 0.7*max(WM_fod_max(:)))); % 20/12/2017
            %             for ij=1:size(WM_fod_normalized,2)
            %                 WM_fod_normalized(:,ij) = WM_fod_normalized(:,ij) / WM_fod_val;% .* fractions(:,1); % 20/12/2017
            %             end

            output.RSS = reshape(RSS,siz(1:3));
            clear RSS;
            output.FOD = reshape(WM_fod,[siz(1:3),size(WM_fod,2)]);
            clear WM_fod;
            %             output.FOD_norm = reshape(WM_fod_normalized,[siz(1:3),size(WM_fod_normalized,2)]);
            %             clear WM_fod_normalized;

            obj.data.img = permute(obj.data.img,[2 1]);
            obj.data.img = reshape(obj.data.img,siz); % 
        end

        function output = OneVoxelDeconv(obj, Stot,LRKernel, HRKernel,min_bval)
            % Actually perform the deconvolution. All parameters must be set
            % before hand using the dedicated functions. Output is a structure
            % with fields FOD and FOD_norm. FOD_norm is a rescaled version of the FOD
            % meant to be compatible with ExploreDTI-based fiber tractography.
            obj.setInnerShellWeighting(obj.inner_shells_weight);
            NC = obj.n_isotropic;
            ANC = obj.n_anisotropic;
            [~,DeconvMethodCode] = MRTTrack.isSupportedMethod(obj.deconv_method); % -1 = failure; 1 = LSQNONNEG; 2 = DW_RegularizedDeconv; 3 = dRL
            if(DeconvMethodCode == -1)
                warning('Unsupported deconvolution method.');
                return;
            end

            TheTol = 1e-3;
            op_e2 = optimset('TolX',1e-2);
            op = optimset();

            fit_data = true(size(Stot));
            % This normalization assumes there is some b=0s/mm2 data. It is not
            % essential to do it as long as fractions is normalized at the end
            NormFactor = mean(Stot(obj.data.bvals<=min_bval));
            output.S0 = NormFactor;
            Stot = Stot/NormFactor;

            Stot = Stot.*obj.shell_weight; % Weight the lower shells

            if(NC > 0 || ANC > 1)
                piso = zeros(NC,1);
                p_old = Inf; % store the fractions at previous iter

                % The following loop will iterate WM-FOD estimation (with mDRL) and
                % fractions guess (with LSQNONNEG)
                for iter=1:50 % 50 = max number of iterations but usually exits way before

                    DS = max(Stot(fit_data) - obj.weighted_LRKernel(fit_data,end-NC+1:end)*piso,0); % Subtract GM and CSF contributions

                    if(DeconvMethodCode == 4)
                        %                         fODFC = mat_dRL(DS, weighted_LRKernel(:,1:end-NC), 200, obj.NN_L, 8);
                        fODFC = ADT_deconv_RLdamp_1D_noEP(DS, obj.weighted_LRKernel(fit_data,1:end-NC), obj.drl_iters, obj.NN_H);
                    elseif(DeconvMethodCode == 3)
                        %                         fODFC = mat_RL(DS, weighted_LRKernel(:,1:end-NC), 200);
                        fODFC = RichardsonLucy(DS, obj.weighted_LRKernel(fit_data,1:end-NC), obj.drl_iters);
                        %                         elseif(DeconvMethodCode == 2)
                        %                             fODFC = DW_RegularizedDeconv(weighted_LRKernel(fit_data,1:end-NC),DS,op_e2,obj.L2LSQ_reg);
                    elseif(DeconvMethodCode == 2 || DeconvMethodCode == 1)
                        fODFC = lsqnonneg(obj.weighted_LRKernel(fit_data,1:end-NC),DS);
                    elseif(DeconvMethodCode == 5)
                        fODFC = ADT_deconv_RLdamp_1D_noEP_white(DS, obj.weighted_LRKernel(fit_data,1:end-NC), obj.drl_iters);
                    elseif(DeconvMethodCode == 6)
                        fODFC = sph_deconv_tv_motor_all(DS, obj.weighted_LRKernel(fit_data,1:end-NC), zeros(size(obj.weighted_LRKernel,2)-NC,1), obj.drl_iters, 1, 'SMF-SENSE-based');
                    elseif(DeconvMethodCode == 7)
                        fODFC = lasso(obj.weighted_LRKernel(fit_data,1:end-NC),DS,'Alpha',1,'Lambda',obj.L2LSQ_reg);
                    end
                    % This line is quite tricky. It comes out of some trial and error
                    % but actually has the intent of eliminating 1) small contributions
                    % 2) flat - spread fods enforcing sparsity.
                    % Without this line the code DOESN'T work. (WM is over-estimated).
                    fODFC(fODFC < median(fODFC)) = 0;
                    % Build a dictionary to fit the complete signal (Stot)
                    Y = [LRKernel(:,1:end-NC)*fODFC LRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                    M = max(Y(:,1));
                    if(M > 0) % i.e. if the FOD is non-zero
                        Y(:,1) = Y(:,1)/M; % Normalize WM signal
                        %                               fODFC = fODFC / sum(fODFC);
                    end

                    p = lsqnonneg(Y(fit_data,:),Stot(fit_data)./obj.shell_weight(fit_data)); % Compute the signal fractions
                    piso = p(end-NC+1:end);

                    % if nothing changed compared to previous iter, exit. (Tol may need to be
                    % adjusted)
                    if(sum(abs(p-p_old) < TheTol) == 3)
                        break
                    end
                    p_old = p;

                    if(obj.robust_deconv == 1 && sum(fit_data == 0) < max_discarded_vols)
                        residuals = Stot./obj.shell_weight - Y*p;
                        fit_data = abs(residuals - median(residuals)) < obj.robust_kappa*1.4826*mad(residuals,1);
                    end

                end

                % NEW FINAL STEP 05/02/2018
                DS = max(Stot(fit_data) - obj.weighted_HRKernel(fit_data,end-NC+1:end)*piso,0); % Subtract GM and CSF contributions

                if(DeconvMethodCode == 4)
                    %                     fODF = mat_dRL(DS, weighted_HRKernel(:,1:end-NC),200, obj.NN_H, 8);
                    fODF = ADT_deconv_RLdamp_1D_noEP(DS, obj.weighted_HRKernel(fit_data,1:end-NC),obj.drl_iters, obj.NN_H);
                elseif(DeconvMethodCode == 3)
                    %                     fODF = mat_RL(DS, weighted_HRKernel(:,1:end-NC), 200);
                    fODF = RichardsonLucy(DS, obj.weighted_HRKernel(fit_data,1:end-NC), obj.drl_iters);
                elseif(DeconvMethodCode == 2)
                    fODF = DW_RegularizedDeconv(obj.weighted_HRKernel(fit_data,1:end-NC),DS,op, obj.L2LSQ_reg);
                elseif(DeconvMethodCode == 1)
                    fODF = lsqnonneg(obj.weighted_HRKernel(fit_data,1:end-NC),DS);
                elseif(DeconvMethodCode == 5)
                    fODF = ADT_deconv_RLdamp_1D_noEP_white(DS, obj.weighted_HRKernel(fit_data,1:end-NC), obj.drl_iters);
                elseif(DeconvMethodCode == 7)
                    fODF = lasso(obj.weighted_HRKernel(fit_data,1:end-NC),DS,'Alpha',1,'Lambda',obj.L2LSQ_reg);
                end
                fODFC = fODF;
                fODFC(fODFC < median(fODFC)) = 0;
                Y = [HRKernel(:,1:end-NC)*fODFC HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                if(obj.n_anisotropic == 1)
                    M = max(Y(:,1));
                    if(M > 0) % i.e. if the FOD is non-zero
                        %                             fODFC = fODFC / sum(fODFC);
                        Y(:,1) = Y(:,1)/M; % Normalize WM signal
                    end
                else
                    Y = [];
                    cindex = 1;
                    for kANC = 1:ANC
                        %         if(sum(Y(:,1)) > 0) % i.e. if the FOD is non-zero
                        % changed for multi FOD
                        NewCol = obj.weighted_HRKernel(:,cindex:cindex+obj.nDirections-1)*fODFC(cindex:cindex+obj.nDirections-1);
                        if(sum(NewCol>0) > 0)
                            NewCol = NewCol/max(NewCol);
                        end
                        Y = [Y NewCol];
                        cindex = cindex+obj.nDirections;
                    end
                    Y = [Y obj.weighted_HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                end
                p = lsqnonneg(Y(fit_data,:),Stot(fit_data)./obj.shell_weight(fit_data),op_e2); % Compute the signal fractions
                output.RSS = sum((Stot-Y*p).^2);

                output.fractions = p/sum(p);
                output.WM_fod = single(fODF);
            else
                NormFactor = mean(Stot(obj.data.bvals<=min_bval));
                Stot = Stot/NormFactor;

                output.fractions = 1;
                FOD = ADT_deconv_RLdamp_1D_noEP(Stot, obj.weighted_HRKernel,obj.drl_iters, obj.NN_H);
                output.WM_fod = single(FOD);
                output.RSS = sum((Stot-obj.weighted_HRKernel*FOD).^2);
            end
        end

        function output = OneVoxelDeconvTest2(obj, Stot,LRKernel, HRKernel,min_bval)
            % Actually perform the deconvolution. All parameters must be set
            % before hand using the dedicated functions. Output is a structure
            % with fields FOD and FOD_norm. FOD_norm is a rescaled version of the FOD
            % meant to be compatible with ExploreDTI-based fiber tractography.
            obj.setInnerShellWeighting(obj.inner_shells_weight);
            NC = obj.n_isotropic;
            ANC = obj.n_anisotropic;
            [~,DeconvMethodCode] = MRTTrack.isSupportedMethod(obj.deconv_method); % -1 = failure; 1 = LSQNONNEG; 2 = DW_RegularizedDeconv; 3 = dRL
            if(DeconvMethodCode == -1)
                warning('Unsupported deconvolution method.');
                return;
            end

            TheTol = 1e-3;
            op_e2 = optimset('TolX',1e-2);
            op = optimset();

            fit_data = true(size(Stot));
            % This normalization assumes there is some b=0s/mm2 data. It is not
            % essential to do it as long as fractions is normalized at the end
            NormFactor = mean(Stot(obj.data.bvals<=min_bval));
            output.S0 = NormFactor;
            Stot = Stot/NormFactor;

            Stot = Stot.*obj.shell_weight; % Weight the lower shells

            K = obj.weighted_LRKernel(fit_data,1:end-NC);
            KT = obj.weighted_LRKernel(fit_data,1:end-NC)';

            if(NC > 0 || ANC > 1)
                piso = zeros(NC,1);
                p_old = Inf; % store the fractions at previous iter

                % The following loop will iterate WM-FOD estimation (with mDRL) and
                % fractions guess (with LSQNONNEG)
                fODF = ones(size(obj.LRKernel,2)-NC,1);
                for iter=1:50 % 50 = max number of iterations but usually exits way before

                    DS = max(Stot(fit_data) - obj.weighted_LRKernel(fit_data,end-NC+1:end)*piso,0); % Subtract GM and CSF contributions

                    if(DeconvMethodCode == 4)
                        %                         fODFC = mat_dRL(DS, weighted_LRKernel(:,1:end-NC), 200, obj.NN_L, 8);
                        fODFC = ADT_deconv_RLdamp_1D_noEP_init(fODF,DS, K, KT, obj.drl_iters, obj.NN_H);
                    elseif(DeconvMethodCode == 3)
                        %                         fODFC = mat_RL(DS, weighted_LRKernel(:,1:end-NC), 200);
                        fODFC = RichardsonLucy(DS, obj.weighted_LRKernel(fit_data,1:end-NC), obj.drl_iters);
                        %                         elseif(DeconvMethodCode == 2)
                        %                             fODFC = DW_RegularizedDeconv(weighted_LRKernel(fit_data,1:end-NC),DS,op_e2,obj.L2LSQ_reg);
                    elseif(DeconvMethodCode == 2 || DeconvMethodCode == 1)
                        fODFC = lsqnonneg(obj.weighted_LRKernel(fit_data,1:end-NC),DS);
                    end
                    % This line is quite tricky. It comes out of some trial and error
                    % but actually has the intent of eliminating 1) small contributions
                    % 2) flat - spread fods enforcing sparsity.
                    % Without this line the code DOESN'T work. (WM is over-estimated).
                    fODFC(fODFC < median(fODFC)) = 0;
%                     fODFC(fODFC < 0.01 | fODFC < 0.1*max(fODFC)) = 0;
                    % Build a dictionary to fit the complete signal (Stot)
                    Y = [LRKernel(:,1:end-NC)*fODFC LRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                    M = max(Y(:,1));
                    if(M > 0) % i.e. if the FOD is non-zero
                        Y(:,1) = Y(:,1)/M; % Normalize WM signal
                        %                               fODFC = fODFC / sum(fODFC);
                    end

                    p = lsqnonneg(Y(fit_data,:),Stot(fit_data)./obj.shell_weight(fit_data)); % Compute the signal fractions
                    piso = p(end-NC+1:end);

                    % if nothing changed compared to previous iter, exit. (Tol may need to be
                    % adjusted)
                    if(sum(abs(p-p_old) < TheTol) == 3)
                        break
                    end
                    p_old = p;

                    if(obj.robust_deconv == 1 && sum(fit_data == 0) < max_discarded_vols)
                        residuals = Stot./obj.shell_weight - Y*p;
                        fit_data = abs(residuals - median(residuals)) < obj.robust_kappa*1.4826*mad(residuals,1);
                    end

                end

                K = obj.weighted_HRKernel(fit_data,1:end-NC);
                KT = obj.weighted_HRKernel(fit_data,1:end-NC)';
                % NEW FINAL STEP 05/02/2018
                DS = max(Stot(fit_data) - obj.weighted_HRKernel(fit_data,end-NC+1:end)*piso,0); % Subtract GM and CSF contributions
                fODF = ones(size(obj.weighted_HRKernel,2)-NC,1);
                if(DeconvMethodCode == 4)
                    fODF = ADT_deconv_RLdamp_1D_noEP_init(fODF, DS, K, KT,obj.drl_iters, obj.NN_H);
                elseif(DeconvMethodCode == 3)
                    fODF = RichardsonLucy(DS, obj.weighted_HRKernel(fit_data,1:end-NC), obj.drl_iters);
                elseif(DeconvMethodCode == 2)
                    fODF = DW_RegularizedDeconv(obj.weighted_HRKernel(fit_data,1:end-NC),DS,op, obj.L2LSQ_reg);
                elseif(DeconvMethodCode == 1)
                    fODF = lsqnonneg(obj.weighted_HRKernel(fit_data,1:end-NC),DS);
                end
                fODFC = fODF;
                fODFC(fODFC < 0.01 | fODFC < 0.1*max(fODFC)) = 0;
                Y = [HRKernel(:,1:end-NC)*fODFC HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                if(obj.n_anisotropic == 1)
                    M = max(Y(:,1));
                    if(M > 0) % i.e. if the FOD is non-zero
                        %                             fODFC = fODFC / sum(fODFC);
                        Y(:,1) = Y(:,1)/M; % Normalize WM signal
                    end
                else
                    Y = [];
                    cindex = 1;
                    for kANC = 1:ANC
                        %         if(sum(Y(:,1)) > 0) % i.e. if the FOD is non-zero
                        % changed for multi FOD
                        NewCol = obj.weighted_HRKernel(:,cindex:cindex+obj.nDirections-1)*fODFC(cindex:cindex+obj.nDirections-1);
                        if(sum(NewCol>0) > 0)
                            NewCol = NewCol/max(NewCol);
                        end
                        Y = [Y NewCol];
                        cindex = cindex+obj.nDirections;
                    end
                    Y = [Y obj.weighted_HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                end
                p = lsqnonneg(Y(fit_data,:),Stot(fit_data)./obj.shell_weight(fit_data),op_e2); % Compute the signal fractions
                output.RSS = sum((Stot-Y*p).^2);

                output.fractions = p/sum(p);
                output.WM_fod = single(fODFC);
            else
                NormFactor = mean(Stot(obj.data.bvals<=min_bval));
                Stot = Stot/NormFactor;

                output.fractions = 1;
                FOD = ADT_deconv_RLdamp_1D_noEP(Stot, obj.weighted_HRKernel,obj.drl_iters, obj.NN_H);
                output.WM_fod = single(FOD);
                output.RSS = sum((Stot-obj.weighted_HRKernel*FOD).^2);
            end
        end

        function output = OneVoxelDeconvTest3(obj, Stot,LRKernel, HRKernel,min_bval)
            % Actually perform the deconvolution. All parameters must be set
            % before hand using the dedicated functions. Output is a structure
            % with fields FOD and FOD_norm. FOD_norm is a rescaled version of the FOD
            % meant to be compatible with ExploreDTI-based fiber tractography.
            obj.setInnerShellWeighting(obj.inner_shells_weight);
            NC = obj.n_isotropic;
            ANC = obj.n_anisotropic;
            [~,DeconvMethodCode] = MRTTrack.isSupportedMethod(obj.deconv_method); % -1 = failure; 1 = LSQNONNEG; 2 = DW_RegularizedDeconv; 3 = dRL
            if(DeconvMethodCode == -1)
                warning('Unsupported deconvolution method.');
                return;
            end

            TheTol = 1e-3;
            op_e2 = optimset('TolX',1e-2);
            op = optimset();

            fit_data = true(size(Stot));
            % This normalization assumes there is some b=0s/mm2 data. It is not
            % essential to do it as long as fractions is normalized at the end
            NormFactor = mean(Stot(obj.data.bvals<=min_bval));
            output.S0 = NormFactor;
            Stot = Stot/NormFactor;

            Stot = Stot.*obj.shell_weight; % Weight the lower shells

            K = obj.weighted_LRKernel(fit_data,1:end-NC);
            KT = obj.weighted_LRKernel(fit_data,1:end-NC)';

            if(NC > 0 || ANC > 1)

                % The following loop will iterate WM-FOD estimation (with mDRL) and
                % fractions guess (with LSQNONNEG)

                K = obj.weighted_HRKernel(fit_data,1:end-NC);
                KT = obj.weighted_HRKernel(fit_data,1:end-NC)';
                % NEW FINAL STEP 05/02/2018
                fODF = ones(size(obj.weighted_HRKernel,2)-NC,1);
                if(DeconvMethodCode == 4)
                    fODF = ADT_deconv_RLdamp_1D_noEP_init(fODF, Stot, K, KT,obj.drl_iters, obj.NN_H);
                elseif(DeconvMethodCode == 3)
                    fODF = RichardsonLucy(Stot, obj.weighted_HRKernel(fit_data,1:end-NC), obj.drl_iters);
                elseif(DeconvMethodCode == 2)
                    fODF = DW_RegularizedDeconv(obj.weighted_HRKernel(fit_data,1:end-NC),Stot,op, obj.L2LSQ_reg);
                elseif(DeconvMethodCode == 1)
                    fODF = lsqnonneg(obj.weighted_HRKernel(fit_data,1:end-NC),Stot);
                end
                fODFC = fODF;
                fODFC(fODFC < 0.01 | fODFC < 0.1*max(fODFC)) = 0;
                Y = [HRKernel(:,1:end-NC)*fODFC HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                if(obj.n_anisotropic == 1)
                    M = max(Y(:,1));
                    if(M > 0) % i.e. if the FOD is non-zero
                        %                             fODFC = fODFC / sum(fODFC);
                        Y(:,1) = Y(:,1)/M; % Normalize WM signal
                    end
                else
                    Y = [];
                    cindex = 1;
                    for kANC = 1:ANC
                        %         if(sum(Y(:,1)) > 0) % i.e. if the FOD is non-zero
                        % changed for multi FOD
                        NewCol = obj.weighted_HRKernel(:,cindex:cindex+obj.nDirections-1)*fODFC(cindex:cindex+obj.nDirections-1);
                        if(sum(NewCol>0) > 0)
                            NewCol = NewCol/max(NewCol);
                        end
                        Y = [Y NewCol];
                        cindex = cindex+obj.nDirections;
                    end
                    Y = [Y obj.weighted_HRKernel(:,end-NC+1:end)]; % 3 columns (WM-SIGNAL GM-SIGNAL CSF-SIGNAL)
                end
                p = lsqnonneg(Y(fit_data,:),Stot(fit_data)./obj.shell_weight(fit_data),op_e2); % Compute the signal fractions
                output.RSS = sum((Stot-Y*p).^2);

                output.fractions = p/sum(p);
                output.WM_fod = single(fODFC);
            else
                NormFactor = mean(Stot(obj.data.bvals<=min_bval));
                Stot = Stot/NormFactor;

                output.fractions = 1;
                FOD = ADT_deconv_RLdamp_1D_noEP(Stot, obj.weighted_HRKernel,obj.drl_iters, obj.NN_H);
                output.WM_fod = single(FOD);
                output.RSS = sum((Stot-obj.weighted_HRKernel*FOD).^2);
            end
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

        function SaveOutputToNii(SpherDec,output,file_prefix)
            % Save the content of a deconvolution data structure to nifti. SpherDec
            % is an instance of this class, output is the structure returned
            % from PerformDeconv. file_prefix is the name without extension of
            % the desired outputs.
            lmax = SpherDec.Lmax_FOD;
            super_scheme = gen_scheme(SpherDec.nDirections,lmax); % the reconstruction scheme. Change 300 to any number
            %             super_scheme.vert = super_scheme.vert(:,[2 1 3]);
            %             super_scheme.vert(:,3) = -super_scheme.vert(:,3);

            sh = SH(lmax,super_scheme.vert);
            max_bv = max(round(SpherDec.data.bvals));
            shell = abs(SpherDec.data.bvals - max_bv) <= 0.01*max_bv;

            if(SpherDec.n_anisotropic == 1)
                % Single FOD case
                test_coef = sh.coef(output.FOD(1,1,1,:));
                fod = zeros([size(output.FOD(:,:,:,1)) length(test_coef)],'single');
                S = sum(SpherDec.HRKernel(shell,1))*SpherDec.nDirections/300;
                output.FOD = output.FOD * S;
                %                 for ij=1:size(F,4)
                %                     F(:,:,:,ij) = F(:,:,:,ij).*output.fractions(:,:,:,1);
                % %                     fod(:,:,:,ij) = fod(:,:,:,ij);
                %                 end
                for sl=1:10:size(output.FOD,3)
                    indices = sl:min(sl+9,size(output.FOD,3));
                    fod(:,:,indices,:) = sh.coef(output.FOD(:,:,indices,:));
                end
                %             FOD_max = max(output.FOD,[],4);
                %             FOD_scaled = output.FOD;
                %             FOD_val = mean(FOD_max(output.fractions(:,:,:,1) > 0.7*max(FOD_max(:)))); % 20/12/2017
                %             for ij=1:size(FOD_scaled,4)
                %                 FOD_scaled(:,:,:,ij) = FOD_scaled(:,:,:,ij) / FOD_val;% .* fractions(:,:,:,1); % 20/12/2017
                %             end


                if(isfield(SpherDec.data,'hdr'))
                    DW_SaveVolumeLikeNii(fod,SpherDec.data,[file_prefix '_CSD_FOD'],0);
                    DW_SaveVolumeLikeNii(output.fractions,SpherDec.data,[file_prefix '_fractions'],0);
                else
                    data_struct.img = fod;
                    data_struct.VD = SpherDec.data.VD;
                    MRTQuant.WriteNifti(data_struct,[file_prefix '_CSD_FOD.nii']);
                    data_struct.img = output.fractions;
                    MRTQuant.WriteNifti(data_struct,[file_prefix '_fractions.nii']);
                end
            else
                % multi-FOD case
                for fod_id=1:SpherDec.n_anisotropic
                    idx = 1+(fod_id-1)*SpherDec.nDirections:fod_id*SpherDec.nDirections;
                    fod = sh.coef(output.FOD(:,:,:,idx));
                    S = sum(SpherDec.HRKernel(shell,idx(1)))*SpherDec.nDirections/300;
                    fod = fod*S;
                    %                     for ij=1:size(fod,4)
                    %                         fod(:,:,:,ij) = fod(:,:,:,ij).*output.fractions(:,:,:,fod_id);
                    %                     end

                    if(isfield(SpherDec.data,'hdr'))
                        DW_SaveVolumeLikeNii(fod,SpherDec.data,[file_prefix '_CSD_FOD_' num2str(fod_id)],0);
                    else
                        data_struct.img = fod;
                        data_struct.VD = SpherDec.data.VD;
                        MRTQuant.WriteNifti(data_struct,[file_prefix '_CSD_FOD_' num2str(fod_id) '.nii']);
                    end
                end
                if(isfield(SpherDec.data,'hdr'))
                    DW_SaveVolumeLikeNii(output.fractions,SpherDec.data,[file_prefix '_fractions'],0);
                else
                    data_struct.img = output.fractions;
                    MRTQuant.WriteNifti(data_struct,[file_prefix '_fractions.nii']);
                end

                json.CallFunction = 'SphericalDeconvolution.PerformDeconv';
                json.Description = 'GRL/mFOD spherical deconvolution';
                json.parameters.n_anisotropic = SpherDec.n_anisotropic;
                json.parameters.n_isotropic = SpherDec.n_isotropic;
                json.parameters.deconv_method = SpherDec.deconv_method;
                json.parameters.inner_shells_weight = SpherDec.inner_shells_weight;
                json.parameters.rf_models = SpherDec.rf_models;
                json.parameters.L2LSQ_reg = SpherDec.L2LSQ_reg;
                json.parameters.nDirections = SpherDec.nDirections;
                json.parameters.NN_L = SpherDec.NN_L;
                json.parameters.NN_H = SpherDec.NN_H;

                NiftiIO_basic.WriteJSONDescription('output',file_prefix,'props',json);

            end
        end

    end

    methods(Static)
        function methods = SupportedMethods()
            % List the supported deconvolution methods
            methods = {'LSQ','L2LSQ','RL','dRL','dRL_W','RUMBA','LASSO','CONSTR'};
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

            intersect_mask = MRTQuant.LoadNifti(fraction_file);
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
                        if(~strcmp(mask_mode,'gm_only'))
                            TP = tract(points2keep_l,:);
                            TL = 0;
                            for tpid=2:size(TP,1)
                                TL = TL + norm(TP(tpid,:)-TP(tpid-1,:));
                            end
                            if(TL < Parameters.Length_range(1))
                                continue
                            end
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

        function [init_lambdas,init_K, init_offset] = Eigenval_IsoK_WM_FromData(data,mask,fit_dki,fit_offset)
            % Estimate the tensor eigenvalues and isotropic kurtosis from the data, in a mask where FA >= 0.7
            data.img = single(data.img);
            [sx,sy,sz,st] = size(data.img);
            signal_stack = reshape(data.img,sx*sy*sz,st);

            if(exist('fit_dki','var') < 1)
                fit_dki = 1;
            end
            if(exist('fit_offset','var') < 1)
                fit_offset = 0;
                init_offset = NaN;
            end

            if(fit_dki == 1)
                Np = 8;
            else
                Np = 7;
            end

            [G,WG,GoodIndexes] = DW_BuildDTMat(data,unique(data.bvals),1);
            Gt = G;
            if(fit_dki == 1)
                G = [Gt(:,1:6) 1/6*(data.bvals*1e-3).^2 Gt(:,end)]; % Extend for isotropic kurtosis
            end
            msk = sum(signal_stack==0,2) == 0;
            p = zeros(length(signal_stack),Np);
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
                p1 = EDTI_Library.E_DTI_WLLS_WW(S,G); % This uses the weighted least squares as implemented in ExploreDTI
                p(ij,:) = p1;%G\log(S);
                if(sum(~isfinite(p(ij,:))) > 0)
                    continue
                end
                autoval = eig(D2Dtensor(p(ij,:)));
                autovals(ij,:) = sort(autoval,'descend');
                K(ij) = 1e-6*p1(end-1)/mean(autoval.^2);
            end

            p = reshape(p,sx,sy,sz,Np);
            DT = EDTI_Library.E_DTI_DWI_mat2cell(p(:,:,:,1:6));
            [~,m_FA,~,~,m_lambdas] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
            m_lambdas = real(m_lambdas);
            K = real(K);
            m_FA = m_FA / sqrt(3);

            U = m_FA(:) > 0.7 & mask(:) > 0;

            if(fit_offset == 0)
                m_lambdas = reshape(m_lambdas,sx*sy*sz,3);
                init_lambdas = median(m_lambdas(U > 0,:));
                init_lambdas(2:3) = mean(init_lambdas(2:3));

                if(fit_dki == 1)
                    init_K = mean(K(U(:) > 0 & K(:) > 0 & K(:) < 4));
                else
                    init_K = 0;
                end

                disp(['Calibrated lambdas are:' num2str(init_lambdas)]);
                disp(['Calibrated K is:' num2str(init_K)]);
            else
                
                pt = reshape(p,sx*sy*sz,Np);
                high_fa_voxels = find(U>0);
                high_fa_voxels = high_fa_voxels(randperm(length(high_fa_voxels)));
                high_fa_voxels = high_fa_voxels(1:min(5000,length(high_fa_voxels)));
                
                est = zeros(length(high_fa_voxels),9);
                params.bvals = data.bvals;
%                 params.bvecs = data.bvecs;
                params.Gt = Gt(:,1:6);
                x0 = [exp(mean(pt(high_fa_voxels,end))) zeros(1,6) 0 0];
                lb = [x0(1)/2 -3e-3*ones(1,6) 0 0];
                ub = [x0(1)*10 3e3*ones(1,6) 3*3e3 x0(1)*10];
                signal_stack = signal_stack(high_fa_voxels,:);
                for x=1:size(est,1)
                   params.S = signal_stack(x,:)';
                   est(x,:) = lsqnonlin(@dki_offset,x0,lb,ub,optimset('display','off'),params); 
                end
                DT = EDTI_Library.E_DTI_DWI_mat2cell(reshape(est(:,2:7),[1 1 size(est,1) 6]));
                [~,m_FA,~,~,m_lambdas] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
                init_lambdas = mean(squeeze(m_lambdas));
                MD = mean(squeeze(m_lambdas),2);
                init_K = mean(est(:,8)./MD.^2);
                init_offset = mean(est(:,9));
                % init_lambdas to be added
            end

            function out = dki_offset(x0,params)
                    s0 = x0(1);
                    D = x0(2:7)';%[x0(2) x0(3) x0(4)
                        %x0(3) x0(5) x0(6)
                        %x0(4) x0(6) x0(7)];
%                     [eigvec,eigval] = eig(D);
%                     MD = mean(diag(eigval));
                    MK = x0(8);
                    offset = x0(9);
%                     Sg = s0*exp(-params.bvals.*diag(params.bvecs*...
%                         D*params.bvecs')+1/6*params.bvals.^2*MD^2*MK)+offset;
                    Sg = s0*(exp(params.Gt*D+1/6*params.bvals.^2*MK)+offset);

                    out = double(params.S - Sg);
            end            
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
                [bvals,bvecs] = MRTQuant.bval_bvec_from_b_Matrix(bmat);
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
            data.bvecs(:,1) = -data.bvecs(:,1);
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
                [bvals,bvecs] = MRTQuant.bval_bvec_from_b_Matrix(bmat);
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
            % nvert: if gradients are not specified - the number of points
            % on the unit sphere (default 300)
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
            nvert = GiveValueForName(varargin,'nvert');
            if(isempty(shcoeffs))
                error('Missing mandatory input shcoeffs');
            end

            VD = [];
            if(ischar(shcoeffs))
                img = MRTQuant.LoadNifti(shcoeffs);
                shcoeffs = img.img;
                VD = img.VD;
                clear img;
            end

            if(isempty(bvals) && isempty(bmat))
                if(isempty(nvert))
                    nvert = 300;
                end
                bvals = 3000*ones(nvert,1);
                Q = gen_scheme(nvert,4);
                bvecs = Q.vert;
            end

            if(isempty(basis_in))
                error('Missing mandatory argument basis_in');
            end
            if(isempty(basis_out))
                error('Missing mandatory argument basis_out');
            end
            if(~isempty(bmat))
                if(ischar(bmat))
                    bmat = load(bmat);
                end
            end

            data = SphericalDeconvolution.SH_2_dMRI('bvals',bvals,...
                'bvecs',bvecs,'basis',basis_in,'bmat',bmat,'shcoeffs',shcoeffs,'lmax',lmax);
            SHCoeffs = SphericalDeconvolution.dMRI_2_SH('data',data,'basis',...
                basis_out,'lmax',lmax);

            if(~isempty(VD) && ~isempty(output))
                out.VD = VD;
                out.img = SHCoeffs;
                clear data SHCoeffs
                MRTQuant.WriteNifti(out,output);
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
                DATA = MRTQuant.LoadNifti(fod_file);
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
                C = MRTQuant.LoadNifti(fod);
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
            parallel = GiveValueForName(varargin,'parallel');
            if(isempty(parallel))
                parallel = 1;
            end

            SHPrecomp.init(lmax,300);
            c = parcluster('local');
            max_cores = c.NumWorkers;
            if(parallel == 0)
                max_cores = 1;
            end

            if(parallel == 1)
                [peak_dir,peak_amp] = SHPrecomp.all_peaks_parallel(fod,peak_threshold,max_peaks,max_cores,lmax,300);
            else
                [peak_dir,peak_amp] = SHPrecomp.all_peaks(fod,peak_threshold,max_peaks,max_cores,lmax,300);
            end

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
            %             shift = round(size(T.TractMask).*T.VDims)/2;
            shift = (size(T.TractMask).*T.VDims)/2;
            for ij=1:length(T.Tracts)
                TR = T.Tracts{ij}(:,[2 1 3]);
                TR(:,1) = -TR(:,1) + shift(2);
                TR(:,2) = -TR(:,2) + shift(1);
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

            ref = MRTQuant.LoadNifti(nii_file);
            VDims = ref.VD;
            TractMask = ones(size(ref.img));
            TractMask = permute(TractMask,[2 1 3]);
            %             shift = (size(TractMask).*VDims)/2;
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

        function ConvertMatTractography2TRK(varargin)
            % Converts MAT tractography 2 the TRK format
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

            TR = load(mat_file,'Tract','VDims','TractMask');

        end

        function ConvertTRKTractography2Mat(varargin)
            % Converts TRK tractography 2 the ExploreDTI format
            % input arguments:
            % trk_file: the input tractography trk file
            % nii_file: a reference .nii file in the same space of the tck
            % output: the output MAT file
            trk_file = GiveValueForName(varargin,'trk_file');
            if(isempty(vtk_file))
                error('Missing mandatory argument trk_file');
            end
            output = GiveValueForName(varargin,'output');
            if(isempty(output))
                error('Missing mandatory argument output');
            end
            nii_file = GiveValueForName(varargin,'nii_file');
            if(isempty(nii_file))
                error('Missing mandatory argument nii_file');
            end

            [header,tracks] = trk_read(trk_file);
            SEG = MRTQuant.LoadNifti(nii_file);

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
                tracks(tid).matrix = tracks(tid).matrix(:,[2 1 3]);
                % These lines are needed for .trk saved from TrackVis, but not for .trk converted
                % from E_DTI - CHECK THE DENSITY MAPS EVERY TIME!!!
                tracks(tid).matrix(:,1) = size(SEG.img,1)-tracks(tid).matrix(:,1);
                tracks(tid).matrix(:,2) = size(SEG.img,2)-tracks(tid).matrix(:,2);
                PL = tracks(tid).matrix;
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
                    %                    cTractFA(point) = DIFF_PROP.FA(PC(1),PC(2),PC(3))/sqrt(3);
                    %                    cTractFE(point,:) = abs(squeeze(DIFF_PROP.FE(PC(1),PC(2),PC(3),:)));
                    %                    cTractLambdas(point,:) = squeeze(DIFF_PROP.eigval(PC(1),PC(2),PC(3),:));
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

        function [mx,my,mz] = ConvertMatTractography2VTK(varargin)
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

        function ConvertVTKTractography2Mat(varargin)
            % Converts VTK tractography 2 the ExploreDTI format
            % input arguments:
            % vtk_file: the input tractography vtk file
            % nii_file: a reference .nii file in the same space of the tck
            % mat_file: a reference .mat file in the same space of the tck
            % (alternative to nii_file)
            % offset: an [x,y,z] vector specificing a 3D translation
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
            mat_file = GiveValueForName(varargin,'mat_file');
            if(isempty(nii_file))
                if(isempty(mat_file))
                    error('Missing mandatory argument nii_file or mat_file');
                else
                    R = load(mat_file,'VDims','FA','FE','eigval');
                    REF.VD = R.VDims;
                    REF.img = R.FA;
                    REF.FA = R.FA;
                    REF.FE = R.FE;
                    REF.eigval = R.eigval;
                    clear R
                end
            else
                REF = MRTQuant.LoadNifti(nii_file);
            end

            offset = GiveValueForName(varargin,'offset');
            if(isempty(offset))
                offset = [0 0 0];
            end

            try
                [Points,Lines] = VTKImport_test(vtk_file);
            catch err
                disp('Error trying to read the VTK file');
                return
            end
            Points(:,1) = Points(:,1) + offset(1);
            Points(:,2) = Points(:,2) + offset(2);
            Points(:,3) = Points(:,3) + offset(3);
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
                if(isfield(REF,'FA'))
                    % These values are retrieved only if specifying
                    % mat_file in place of nii_file
                    for point=1:size(PL,1)
                        P = PL(point,:)./REF.VD;
                        PC = floor(P)+1;
                        cTractFA(point) = REF.FA(PC(1),PC(2),PC(3))/sqrt(3);
                        cTractFE(point,:) = abs(squeeze(REF.FE(PC(1),PC(2),PC(3),:)));
                        cTractLambdas(point,:) = squeeze(REF.eigval(PC(1),PC(2),PC(3),:));
                        cTractMD(point) = mean(cTractLambdas(point,:));
                    end
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

        function FilterTractsANDGate(varargin)
            % filt_type 0: keep if pass once in params.mask1
            % filt_type 1: keep if pass at least twice in params.mask1
            % filt_type 2: keep if start in params.mask1 and end in params.mask2
            % filt_type 3: discard if passing in params.mask1
            % filt_type 4: keep if starts and/or ends in params.mask1
            coptions = varargin;
            base_file = GiveValueForName(coptions,'mat_file');
            if(isempty(base_file))
                error('Need to specify the input .mat file');
            end
            tract_file = GiveValueForName(coptions,'tract_file');
            if(isempty(tract_file))
                error('Need to specify the input Tract.mat file');
            end
            out_file = GiveValueForName(coptions,'output');
            if(isempty(out_file))
                error('Need to specify the output file');
            end
            params = GiveValueForName(coptions,'params');
            if(isempty(out_file))
                error('Need to specify input parameters');
            end
            filt_type = GiveValueForName(coptions,'filt_type');
            if(isempty(filt_type))
                error('Need to specify input parameters');
            end

            load(tract_file);
            load(base_file,'FA','VDims');

            [sx,sy,sz] = size(FA);

            good_tracts = true(length(Tracts),1);

            if(filt_type < 2)
                params.mask2 = params.mask1;
            end

            if(exist('min_length','var') < 1)
                min_length = 2;
            end

            try

                for tid=1:length(Tracts)
                    tract = Tracts{tid};
                    passed_1 = 0;
                    passed_2 = 0;
                    start_1 = 0;
                    start_2 = 0;
                    end_1 = 0;
                    end_2 = 0;

                    for l=1:size(tract,1)
                        point = round(tract(l,:)./VDims);
                        if(params.mask1(point(1),point(2),point(3)))
                            passed_1 = passed_1 + 1;
                            if(l == 1)
                                start_1 = 1;
                            elseif(l == size(tract,1))
                                end_1 = 1;
                            end
                        end
                        if(isfield(params,'mask2'))
                            if(params.mask2(point(1),point(2),point(3)))
                                passed_2 = passed_2 + 1;
                                if(l == 1)
                                    start_2 = 1;
                                elseif(l == size(tract,1))
                                    end_2 = 1;
                                end
                            end
                        end
                    end
                    if(filt_type == 0)
                        if(passed_1 < 1)
                            good_tracts(tid) = false;
                        end
                    elseif(filt_type == 1)
                        if(passed_1 < 2)
                            good_tracts(tid) = false;
                        end
                    elseif(filt_type == 2)
                        if((start_1 < 1 && start_2 < 1) || (end_1 < 1 && ...
                                end_2 < 1))
                            good_tracts(tid) = false;
                        else
                            if((start_1 > 0 && end_2 < 1) || ...
                                    (start_2 > 0 && end_1 < 1))
                                good_tracts(tid) = false;
                            end
                        end
                    elseif(filt_type == 3)
                        if(passed_1 > 0 || passed_2 > 0)
                            good_tracts(tid) = false;
                        end
                    elseif(filt_type == 4)
                        if(start_1 < 1 && end_1 < 1)
                           good_tracts(tid) = false; 
                        end
                    end

                    if(TractL{tid} < min_length)
                        good_tracts(tid) = false;
                    end

                end

                Tracts = Tracts(good_tracts);
                if(exist('TractsFOD','var') > 0)
                    TractsFOD = TractsFOD(good_tracts);
                else
                    TractsFOD = [];
                end
                TractL = TractL(good_tracts);
                TractFA = TractFA(good_tracts);
                TractFE = TractFE(good_tracts);
                TractAng = TractAng(good_tracts);
                TractGEO = TractGEO(good_tracts);
                TractLambdas = TractLambdas(good_tracts);
                TractMD = TractMD(good_tracts);
                FList = 1:length(Tracts);

                disp(['Gated ' num2str(sum(good_tracts == true)) ' out of ' num2str(length(good_tracts)) ' ' ...
                    sprintf('%.2f',single(sum(good_tracts == true))/single(length(good_tracts))) '%']);

                save(out_file,'Tracts','TractsFOD','TractL','TractFA','TractFE','TractFE',...
                    'TractAng','TractGEO','TractLambdas','TractMD','TractMask','FList','VDims');

            catch err
                disp(err.message);
            end
        end

        function ConcatenateMATTractographies(varargin)
            % filt_type 0: pass once in params.mask1
            % filt_type 1: pass at least twice in params.mask1
            % filt_type 2: start in params.mask1 and end in params.mask2
            % filt_type 3: discard if passing in params.mask1
            coptions = varargin;
            mat_file1 = GiveValueForName(coptions,'mat_file1');
            if(isempty(mat_file1))
                error('Need to specify at least 2 tractography mat files .mat file');
            end
            output = GiveValueForName(coptions,'output');
            if(isempty(output))
                error('Need to specify the output .mat file');
            end
            TB = load(mat_file1);
            idx = 2;
            while true
                mat_file1 = GiveValueForName(coptions,['mat_file' num2str(idx)]);
                if(isempty(mat_file1))
                    if(idx == 2)
                        error('Need to specify at least 2 tractography mat files .mat file');
                    else
                        break
                    end
                else
                    TA = load(mat_file1);
                    if(any(TA.VDims ~= TB.VDims))
                        error('Incompatible voxel sizes');
                    end
                    NL = length(TA.Tracts);
                    TB.TractAng(end+1:end+NL) = TA.TractAng;
                    TB.TractFA(end+1:end+NL) = TA.TractAng;
                    TB.TractFE(end+1:end+NL) = TA.TractAng;
                    TB.TractGEO(end+1:end+NL) = TA.TractAng;
                    TB.TractL(end+1:end+NL) = TA.TractAng;
                    TB.TractLambdas(end+1:end+NL) = TA.TractAng;
                    TB.TractMD(end+1:end+NL) = TA.TractAng;
                    TB.Tracts(end+1:end+NL) = TA.TractAng;
                    try
                        TB.TractsFOD(end+1:end+NL) = TA.TractAng;
                    catch
                    end
                    TB.FList = 1:length(TB.Tracts);
                end
                idx = idx + 1;
            end
            save(output,'-struct','TB');

        end

        function WMATractographyClustering(varargin)
            % Runs the WMA automatic clustering for fiber tractography
            % see https://github.com/SlicerDMRI/whitematteranalysis for
            % more details. PMID (29920375)
            % Input arguments:
            % mat_file: the reference .MAT file with the DTI/DKI fit (.mat)
            % tract_file: the input tractography file in ExploreDTI (.mat)
            % format
            % output: the output folder containing the clustered
            % tractogaphy
            % keep_temp: whether to keep intermediate temporary files
            % (default false)

            global MRIToolkit;

            mat_file = GiveValueForName(varargin,'mat_file');
            if(isempty(mat_file))
                error('Missing mandatory argument mat_file');
            end

            tract_file = GiveValueForName(varargin,'tract_file');
            if(isempty(tract_file))
                error('Missing mandatory argument tract_file');
            end
            TractsName = tract_file(1:end-4);

            WMA_Folder = GiveValueForName(varargin,'output');
            if(isempty(WMA_Folder))
                error('Missing mandatory argument output');
            end

            keep_temp = GiveValueForName(varargin,'keep_temp');
            if(isempty(keep_temp))
                keep_temp = false;
            end
            
            SlicerPath = MRIToolkit.SlicerPath;
            if(isempty(SlicerPath) || exist(SlicerPath,'dir') < 1)
                error('Please specifiy a valid path to a 3DSlicer installation for your platform in MRIToolkitDefineLocalVars.m');
            end
            SlicerLocation = fullfile(SlicerPath,'Slicer');
            if(ispc)
                SlicerLocation = [SlicerLocation '.exe'];
            end
            if(exist(SlicerLocation,'file') < 1)
                error('Please specifiy a valid path to a 3DSlicer installation for your platform in MRIToolkitDefineLocalVars.m');
            end
            AtlasFolder = MRIToolkit.WMAAtlasFolder;
            if(exist(AtlasFolder,'dir') < 1)
                error('Please specifiy a valid path to the WMA atlas, which can be downloaded from https://github.com/SlicerDMRI/whitematteranalysis/blob/master/doc/subject-specific-tractography-parcellation.md');
            end

            if(exist(WMA_Folder,'dir') < 1)
                mkdir(WMA_Folder);
            end

%             [mx,my,mz] = E_DTI_Convert_tracts_mat_2_vtk_lines(tract_file, strrep(tract_file,'.mat','.vtk'));
            [mx,my,mz] = E_DTI_Convert_tracts_mat_2_vtk_lines(tract_file, 'tracts.vtk');

            conda_path = MRIToolkit.Miniconda3Path;
            conda_env = MRIToolkit.Miniconda3Env;
            if(isempty(conda_env))
                conda_env = conda_path;
            end
            if(isempty(conda_path))
                error('Please specify a valid path to a Conda python3 installation in MRIToolkitDefineLocalVars.m');
            end
            if(ispc)
                conda_path = strrep(conda_path,'\','/');
                base_cmd = ['%windir%\System32\cmd.exe /c "' conda_path '/Scripts/activate.bat ' conda_env ' && '];
            else
                % TODO: account for ZSH on MacOS
                [~,shell] = system('echo $0');
                if(contains(shell,'bash'))
                    base_cmd = ['source ~/.bashrc; conda activate ' conda_env];
                elseif(contains(shell,'zsh'))
                    base_cmd = ['source ~/.zshrc; conda activate ' conda_env];
                else
                    error('Unsupported system shell');
                end
            end

            % cmd = [base_cmd ';wm_quality_control_tract_overlap.py ' ...
            %     AtlasFolder '/ORG-800FC-100HCP/atlas.vtp ' ...
            %     TractsName '.vtk ' WMA_Folder '/QC/InputTractOverlap/'];
            % system(cmd);
            temp_folders = {};

            if(exist([WMA_Folder '/TractRegistration/'],'dir') < 1)
                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_register_to_atlas_new.py']);
                    cmd = [base_cmd ' python ' cmd_loc ...
                        ' -mode rigid_affine_fast ' ...
                        'tracts.vtk ' ...
                        [AtlasFolder '/ORG-RegAtlas-100HCP/registration_atlas.vtk '] ...
                        [WMA_Folder '/TractRegistration']];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';wm_register_to_atlas_new.py ' ...
                        ' -mode rigid_affine_fast ' ...
                        'tracts.vtk ' ...
                        AtlasFolder '/ORG-RegAtlas-100HCP/registration_atlas.vtk ' ...
                        WMA_Folder '/TractRegistration/'];
                end
                system(cmd);
            end

            if(exist([WMA_Folder '/AnatomicalTracts/'],'file') < 1)

                reg_file = dir(fullfile(WMA_Folder,'TractRegistration','*','output_tractography','*.vtk'));
                prop_name = reg_file.name;
                reg_file = fullfile(reg_file.folder,reg_file.name);

                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_cluster_from_atlas.py']);
                    cmd = [base_cmd ' python ' cmd_loc ...
                        reg_file ' '  ...
                        AtlasFolder '/ORG-800FC-100HCP/ ' ...
                        WMA_Folder '/FiberClustering/InitialClusters'];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';wm_cluster_from_atlas.py ' ...
                        reg_file ' '  ...
                        AtlasFolder '/ORG-800FC-100HCP/ ' ...
                        WMA_Folder '/FiberClustering/InitialClusters'];
                end
                system(cmd);
                temp_folders(end+1) = {[WMA_Folder '/FiberClustering']};

                cluster_folder = dir(fullfile(WMA_Folder,'FiberClustering','InitialClusters','*reg'));
                cluster_folder = fullfile(cluster_folder.folder,cluster_folder.name);

                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_cluster_remove_outliers.py']);
                    cmd = [base_cmd ' python ' cmd_loc ...
                        cluster_folder ' ' ...
                        AtlasFolder '/ORG-800FC-100HCP/ ' ...
                        WMA_Folder '/FiberClustering/OutlierRemovedClusters/'];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';wm_cluster_remove_outliers.py ' ...
                        cluster_folder ' ' ...
                        AtlasFolder '/ORG-800FC-100HCP/ ' ...
                        WMA_Folder '/FiberClustering/OutlierRemovedClusters/'];
                end
                system(cmd);

                cluster_folder = dir(fullfile(WMA_Folder,'FiberClustering','OutlierRemovedClusters','*reg*'));
                cluster_folder = fullfile(cluster_folder.folder,cluster_folder.name);

                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_assess_cluster_location_by_hemisphere.py']);
                    cmd = [base_cmd ' python ' cmd_loc ...
                        cluster_folder ' -clusterLocationFile ' ...
                        AtlasFolder '/ORG-800FC-100HCP/cluster_hemisphere_location.txt '];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';wm_assess_cluster_location_by_hemisphere.py ' ...
                        cluster_folder ' -clusterLocationFile ' ...
                        AtlasFolder '/ORG-800FC-100HCP/cluster_hemisphere_location.txt '];
                end
                system(cmd);

                tfm_file = dir(fullfile(WMA_Folder,'TractRegistration','*','output_tractography','*.tfm'));
                tfm_file = fullfile(tfm_file.folder,tfm_file.name);

                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_harden_transform.py']);
                    cmd = [base_cmd ' python ' cmd_loc ' -i -t ' ...
                        tfm_file ' ' cluster_folder '/ '...
                        WMA_Folder '/FiberClustering/TransformedClusters/clustered_data/ ' ...
                        SlicerLocation '"'];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:' SlicerPath ';wm_harden_transform.py -i -t ' ...
                        tfm_file ' ' cluster_folder '/ '...
                        WMA_Folder '/FiberClustering/TransformedClusters/clustered_data/ ' ...
                        SlicerLocation];
                end
                system(cmd);

                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_separate_clusters_by_hemisphere.py']);
                    cmd = [base_cmd ' python ' cmd_loc ...
                        WMA_Folder '/FiberClustering/TransformedClusters/clustered_data/ ' ...
                        WMA_Folder '/FiberClustering/SeparatedClusters/ "'];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';wm_separate_clusters_by_hemisphere.py ' ...
                        WMA_Folder '/FiberClustering/TransformedClusters/clustered_data/ ' ...
                        WMA_Folder '/FiberClustering/SeparatedClusters/ '];
                end
                system(cmd);

                if(ispc)
                    [a,cmd_loc] = system([base_cmd ' where wm_append_clusters_to_anatomical_tracts.py']);
                    cmd = [base_cmd ' python ' cmd_loc ...
                        WMA_Folder '/FiberClustering/SeparatedClusters/ ',...
                        AtlasFolder '/ORG-800FC-100HCP/ ' ...
                        WMA_Folder '/AnatomicalTracts/'];
                    cmd = strrep(cmd,newline,' ');
                else
                    cmd = [base_cmd ';wm_append_clusters_to_anatomical_tracts.py ' ...
                        WMA_Folder '/FiberClustering/SeparatedClusters/ ',...
                        AtlasFolder '/ORG-800FC-100HCP/ ' ...
                        WMA_Folder '/AnatomicalTracts/'];
                end
                system(cmd);

            end

            mkdir(fullfile(WMA_Folder,'AnatomicalTracts_VTK'));
            temp_folders(end+1) = {[WMA_Folder '/AnatomicalTracts_VTK']};

            if(ispc)
                [a,cmd_loc] = system([base_cmd ' where wm_vtp2vtk.py']);
                cmd = [base_cmd ' python ' cmd_loc ...
                    WMA_Folder '/AnatomicalTracts/ ' ...
                    WMA_Folder '/AnatomicalTracts_VTK'];
                cmd = strrep(cmd,newline,' ');
            else
                cmd = [base_cmd ';wm_vtp2vtk.py ' ...
                    WMA_Folder '/AnatomicalTracts/ ' ...
                    WMA_Folder '/AnatomicalTracts_VTK'];
            end
            system(cmd);
            temp_folders(end+1) = {[WMA_Folder '/AnatomicalTracts']};

            vtk_files = dir(fullfile(WMA_Folder,'AnatomicalTracts_VTK','_vtk','*.vtk'));
            output_directory = fullfile(WMA_Folder,'AnatomicalTracts_EDTI');
%             output_directory_density = fullfile(WMA_Folder,'AnatomicalTracts_EDTI_density');
            if(exist(output_directory,'dir') < 1)
                mkdir(output_directory);
            end
%             if(exist(output_directory_density,'dir') < 1)
%                 mkdir(output_directory_density);
%             end

            LM = load([TractsName '.mat'],'TractMask','VDims');
            DIFF_PROP = load(mat_file,'FA','eigval','FE');
            for vid = 1:length(vtk_files)
                try
                    MRTTrack.ConvertVTKTractography2Mat('mat_file',mat_file,...
                        'vtk_file',fullfile(vtk_files(vid).folder,vtk_files(vid).name),...
                        'output',fullfile(output_directory,[vtk_files(vid).name(1:end-4) '.mat']),...
                        'offset',[mx,my,mz]);
%                     E_DTI_DensityMaps(mat_file,fullfile(output_directory,[vtk_files(vid).name(1:end-4) '.mat']),...
%                         fullfile(output_directory_density,[vtk_files(vid).name(1:end-4) '.nii']));
                catch err
                    warning(['Error processing ' vtk_files(vid).name ', continuing']);
                    continue
                end
                %     Neuro.CAT12ApplyDeformation('nii_file',fullfile(output_directory_density,[vtk_files(vid).name(1:end-4) '.nii']),...
                %         'field_file',fullfile(subj_dest_folder,'anat','T1_FP_CAT12','mri','y_T1_FP.nii'));
            end
            if(keep_temp == false)
                for fid=1:length(temp_folders)
                    rmdir(temp_folders{fid},'s');
%                     delete(strrep(tract_file,'.mat','.vtk'));
                    delete('tracts.vtk');
                end
            end
        end

        function PerformDTIBased_FiberTracking(varargin)
            % This function performs whole volume deterministic fiber
            % tractography using the first eigenvector of the DTI fit. Possible
            % input arguments are:
            % mat_file: The ExpoloreDTI-like .mat file
            % output: The output tracts (must be .mat)
            %  SeedPointRes: The seeding resolution in mm, as [2 2 2] (default)
            %  StepSize: the step size in mm, as 1 (default)
            %  FAThresh: The FA thredshold to stop tracking, as 0.2000 	(default)
            %  AngleThresh: The angle threshold to stop tracking, as 30 (default)
            %  FiberLengthRange: The mininum - maximum allowed fiber length in mm: [30 500]
            if(isempty(varargin))
                my_help('MRTTrack.PerformDTIBased_FiberTracking');
                return;
            end

            json.CallFunction = 'MRTTrack.PerformDTIBased_FiberTracking';
            json.Description = my_help('MRTTrack.PerformDTIBased_FiberTracking');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end

            json.ReferenceFile = file_in;
            json.ProcessingType = 'FiberTractography';

            filename_out = GiveValueForName(coptions,'output');
            if(isempty(filename_out))
                error('Need to specify the output .mat file');
            end

            option = GiveValueForName(coptions,'SeedPointRes');
            if(isempty(option))
                parameters.SeedPointRes = [2 2 2];
            else
                parameters.SeedPointRes = option;
            end

            option = GiveValueForName(coptions,'StepSize');
            if(isempty(option))
                parameters.StepSize = 1;
            else
                parameters.StepSize = option;
            end

            option = GiveValueForName(coptions,'FAThresh');
            if(isempty(option))
                parameters.FAThresh = 0.2;
            else
                parameters.FAThresh = option;
            end

            option = GiveValueForName(coptions,'AngleThresh');
            if(isempty(option))
                parameters.AngleThresh = 30;
            else
                parameters.AngleThresh = option;
            end

            option = GiveValueForName(coptions,'FiberLengthRange');
            if(isempty(option))
                parameters.FiberLengthRange = [30 500];
            else
                parameters.FiberLengthRange = option;
            end

            json.parameters = parameters;

            EDTI_Library.WholeBrainTrackingDTI_fast(file_in, filename_out, parameters);
            NiftiIO_basic.WriteJSONDescription('output',filename_out(1:end-4),'props',json);
        end

        function CSD_FOD = PerformCSD(varargin)
            % This function compute the fiber orientation distribution (FOD)
            % using constrained spherical deconvolution (CSD) with recursive
            % response function calibration. Input arguments:
            % mat_file: The .mat file of the data in ExploreDTI-like format
            % mscsd: Use multi-shell multi-tissue CSD in place of classic CSD.
            % rc_mask_file: This is an optional argument to enforce the
            % selection of the response function within the provided .nii mask
            % output: The output .nii where the FOD will be stored
            % Lmax: order of the spherical harmonics. Default is 8.
            % T1seg: pveseg file output of FSL FAST containing 3 tissue classes
            %   (WM, GM, CSF)
            % save_sh: 0 or 1. Save the SH coefficients of the data (csd only)
            % rf_dti: use a DTI based response function in place of recursive calibration.
            %         input as "FA MDx10^3(mm2/s)" of the desired response function, e.g "0.7 1.0"
            %         valid only in combination with csd
            if(isempty(varargin))
                my_help('MRTTrack.PerformCSD');
                return;
            end

            json.CallFunction = 'MRTTrack.PerformCSD';
            json.Description = my_help('MRTTrack.PerformCSD');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end

            json.ReferenceFile = file_in;
            json.ProcessingType = 'Quantification';

            rc_mask_file = GiveValueForName(coptions,'rc_mask_file');
            if(isempty(rc_mask_file))
                rc_mask_file = '';
            end

            json.rc_mask_file = rc_mask_file;

            filename_out = GiveValueForName(coptions,'output');
            if(isempty(filename_out))
                filename_out = [file_in(1:end-4) '_CSD_FOD.nii'];
            end

            option = GiveValueForName(coptions,'Lmax');
            if(~isempty(option))
                Lmax = option;
            else
                Lmax = 8;
            end

            json.Lmax = Lmax;

            option = GiveValueForName(coptions,'T1seg');
            if(~isempty(option))
                t1_seg = option;
            else
                t1_seg = '';
            end

            json.t1_seg = t1_seg;

            save_sh = GiveValueForName(coptions,'save_sh');
            if(isempty(save_sh))
                save_sh = 0;
            end

            json.save_sh = save_sh;

            dti_rf = GiveValueForName(coptions,'rf_dti');
            json.dti_rf = dti_rf;

            mscsd = GiveValueForName(coptions,'mscsd');
            if(isempty(mscsd) || mscsd == 0)
                if(isempty(dti_rf))
                    CSD_FOD = EDTI_Library.E_DTI_Get_HARDI_CSD_FOD_RC(file_in,Lmax,rc_mask_file,filename_out,save_sh);
                else
                    pieces = strsplit(dti_rf);
                    sim_fa = str2double(pieces{1});
                    sim_adc = str2double(pieces{2})*1e-3;
                    CSD_FOD = EDTI_Library.E_DTI_Get_HARDI_CSD_FOD(file_in, 1, Lmax,sim_adc, sim_fa,filename_out,save_sh);
                end
            else
                disp('Running MS-CSD (beta mode)...')
                CSD_FOD = EDTI_Library.E_DTI_Get_HARDI_CSD_FOD_RC_MuSh(file_in,Lmax,rc_mask_file,filename_out,t1_seg);
            end
            json.mscsd = mscsd;

            NiftiIO_basic.WriteJSONDescription('output',filename_out(1:end-4),'props',json);
            if(nargout == 0)
                CSD_FOD = [];
            end

        end

        function PerformFODBased_FiberTracking(varargin)
            % This function performs whole volume deterministic fiber
            % tractography of an FOD in spherical harmonics. Possible
            % input arguments are:
            % mat_file: The ExpoloreDTI-like .mat file (for reference)
            % fod_file: The ExpoloreDTI-like FOD .nii file (in SH basis)
            % output: The output tracts (must be .mat)
            % SeedPointRes: The seeding resolution in mm, as [2 2 2] (default)
            % StepSize: the step size in mm, as 1 (default)
            % FODThresh: The FOD thredshold to stop tracking, as 0.1000 	(default)
            % AngleThresh: The angle threshold to stop tracking, as 30 (default)
            % FiberLengthRange: The mininum - maximum allowed fiber length in mm: [30 500]
            % SeedMask: A mask to perform the seeding. If empty, the whole
            %   volume is used
            % Default parameters:
            %  SeedPointRes: [2 2 2]
            %             StepSize: 1
            %          AngleThresh: 30
            %     FiberLengthRange: [30 500]
            %     FODThresh: 0.1

            if(isempty(varargin))
                my_help('MRTTrack.PerformFODBased_FiberTracking');
                return;
            end

            json.CallFunction = 'MRTTrack.PerformFODBased_FiberTracking';
            json.Description = my_help('MRTTrack.PerformFODBased_FiberTracking');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end

            json.ReferenceFile = file_in;
            json.ProcessingType = 'FiberTractography';

            fod_file = GiveValueForName(coptions,'fod_file');
            if(isempty(fod_file))
                error('Need to specify the input FOD file/variable');
            end
            json.fod_file = fod_file;

            filename_out = GiveValueForName(coptions,'output');
            if(isempty(filename_out))
                error('Need to specify the output .mat file');
            end

            option = GiveValueForName(coptions,'SeedPointRes');
            if(isempty(option))
                parameters.SeedPointRes = [2 2 2];
            else
                parameters.SeedPointRes = option;
            end

            option = GiveValueForName(coptions,'StepSize');
            if(isempty(option))
                parameters.StepSize = 1;
            else
                parameters.StepSize = option;
            end

            option = GiveValueForName(coptions,'FODThresh');
            if(isempty(option))
                parameters.blob_T = 0.1;
            elseif(isnumeric(option))
                parameters.blob_T = option;
            elseif(ischar(option) && strcmpi(option,'auto'))
                % to be adjusted
                Npoints = 5000;

                fod = MRTQuant.LoadNifti(fod_file);
                [sx,sy,sz,st] = size(fod.img);

                Lmax = SH.n2lmax(st);
                SHPrecomp.init(Lmax,300);
                fodV = reshape(fod.img,sx*sy*sz,st);
                % gp = find(csf>0);
                gp = find(fod.img(:,:,:,1) > 0);
                P2R = gp(randperm(length(gp)));
                P2R = P2R(1:Npoints);

                fodV = fodV(P2R,:)';

                [peaks,vals] = SHPrecomp.all_peaks(fodV,0.0001,5);

                gvals = cellfun(@length,vals);
                vals = vals(gvals>0);
                vals = cellfun(@(x)x(1),vals);

                FODThresh = mean(vals)*0.2;
                disp(['FOD thresh is ' num2str(FODThresh)]);
                parameters.blob_T = FODThresh;
            end

            option = GiveValueForName(coptions,'AngleThresh');
            if(isempty(option))
                parameters.AngleThresh = 30;
            else
                parameters.AngleThresh = option;
            end

            option = GiveValueForName(coptions,'FiberLengthRange');
            if(isempty(option))
                parameters.FiberLengthRange = [30 500];
            else
                parameters.FiberLengthRange = option;
            end

            option = GiveValueForName(coptions,'SeedMask');
            if(isempty(option))
                parameters.SeedMask = [];
            else
                parameters.SeedMask = option;
            end

            parameters.randp = 0;

            json.TrackingParameters = parameters;

            EDTI_Library.WholeBrainFODTractography(file_in,fod_file,parameters,filename_out);
            NiftiIO_basic.WriteJSONDescription('output',filename_out(1:end-4),'props',json);
        end

        function Perform_mFOD_FiberTracking(varargin)
            % This function performs whole volume deterministic fiber
            % tractography of multiple FODs in spherical harmonics. Possible
            % input arguments are:
            % mat_file: The ExpoloreDTI-like .mat file (for reference)
            % fod_basename: The name prefix used for the mFOD output (see the SphericalDeconvolution class) (in SH basis)
            % output: The output tracts (must be .mat)
            % SeedPointRes: The seeding resolution in mm, as [2 2 2] (default)
            % StepSize: the step size in mm, as 1 (default)
            % FODThresh: The FOD thredshold to stop tracking, as 0.1000 	(default)
            % AngleThresh: The angle threshold to stop tracking, as 30 (default)
            % FiberLengthRange: The mininum - maximum allowed fiber length in mm: [30 500]
            % SeedMask: A mask to perform the seeding. If empty, the whole
            %   volume is used
            % InterpolationMode: how to merge the mulitple FODs. 'linear'
            %   corresponds to mFOD-WS (each FOD weighted by its fraction and
            %   summed), whereas 'majority' tracks the locally larger
            %   FOD

            if(isempty(varargin))
                my_help('MRTTrack.Perform_mFOD_FiberTracking');
                return;
            end

            json.CallFunction = 'MRTTrack.Perform_mFOD_FiberTracking';
            json.Description = my_help('MRTTrack.PerformFODBased_FiberTracking');

            coptions = varargin;
            file_in = GiveValueForName(coptions,'mat_file');
            if(isempty(file_in))
                error('Need to specify the input .mat file');
            end
            json.ReferenceFile = file_in;
            json.ProcessingType = 'FiberTractography';

            fod_basename = GiveValueForName(coptions,'fod_basename');
            if(isempty(fod_basename))
                error('Need to specify the input FOD basename');
            end
            json.fod_basename = fod_basename;

            filename_out = GiveValueForName(coptions,'output');
            if(isempty(filename_out))
                error('Need to specify the output .mat file');
            end

            option = GiveValueForName(coptions,'SeedPointRes');
            if(isempty(option))
                parameters.SeedPointRes = [2 2 2];
            else
                parameters.SeedPointRes = option;
            end

            option = GiveValueForName(coptions,'StepSize');
            if(isempty(option))
                parameters.StepSize = 1;
            else
                parameters.StepSize = option;
            end

            option = GiveValueForName(coptions,'FODThresh');
            if(isempty(option))
                parameters.blob_T = 0.1;
            else
                parameters.blob_T = option;
            end

            option = GiveValueForName(coptions,'AngleThresh');
            if(isempty(option))
                parameters.AngleThresh = 30;
            else
                parameters.AngleThresh = option;
            end

            option = GiveValueForName(coptions,'FiberLengthRange');
            if(isempty(option))
                parameters.FiberLengthRange = [30 500];
            else
                parameters.FiberLengthRange = option;
            end

            option = GiveValueForName(coptions,'SeedMask');
            if(isempty(option))
                parameters.SeedMask = [];
            else
                parameters.SeedMask = option;
            end

            option = GiveValueForName(coptions,'InterpolationMode');
            if(isempty(option))
                interpolation_mode = 'linear';
            else
                interpolation_mode = option;
            end

            parameters.randp = 0;
            json.TrackingParameters = parameters;

            EDTI_Library.WholeBrainTracking_mDRL_fast_exe(file_in, fod_basename, filename_out, parameters, interpolation_mode);

            NiftiIO_basic.WriteJSONDescription('output',filename_out(1:end-4),'props',json);
        end

        function S = SimulateCrossingFibersSignalWithDKI(varargin)
            % This function generates a synthetic signal using a simplified
            % DKI model (DTI + MK only) of K crossing fibers
            % input arguments are:
            % bvals: the b-values (diffusion weightings vector Nx1)
            % bvecs: the b-vecs (gradients matrix Nx3)
            % fractions: the mixture fractions vector (K elements)
            % eigenvalues: Kx3 matrix of the eigenvalues of each fiber
            % mks: the mean kurtosis of each fiber (K elements)
            % angles: the polar angles of each fiber (Kx2 elements)
            % SNR: The SNR to add Rician Noise (Inf = No Noise)

            coptions = varargin;
            bvals = GiveValueForName(coptions,'bvals');
            if(isempty(bvals))
                error('Need to specify bvals');
            end
            bvecs = GiveValueForName(coptions,'bvecs');
            if(isempty(bvecs))
                error('Need to specify bvecs');
            end
            fractions = GiveValueForName(coptions,'fractions');
            if(isempty(fractions))
                error('Need to specify fractions');
            end
            eigenvalues = GiveValueForName(coptions,'eigenvalues');
            if(isempty(eigenvalues))
                error('Need to specify eigenvalues');
            end
            mks = GiveValueForName(coptions,'mks');
            if(isempty(mks))
                error('Need to specify mks');
            end
            angles = GiveValueForName(coptions,'angles');
            if(isempty(angles))
                error('Need to specify angles');
            end
            SNR = GiveValueForName(coptions,'SNR');
            if(isempty(SNR))
                error('Need to specify SNR');
            end

            fractions = fractions / sum(fractions);
            S = 0;
            for ij = 1:length(fractions)
                %                S = S + fractions(ij)/sum(fractions)*MRT_Library.create_signal_multi_tensor_dki(angles(ij,:), 1, eigenvalues(ij,:), bvals, bvecs, 1, Inf, 0, mks(ij));
                S = S + fractions(ij)*EDTI_Library.SimulateSignalDTIIsotropicK(bvals,bvecs,eigenvalues(ij,:),mks(ij),angles(ij,:));
            end

            if(~isinf(SNR))
                S = MRT_Library.AddRicianNoise(S,1/SNR);
            end

        end

        % TO DO: harmonize inputs and write help
        function SubsampleTractsUniform(tract_file,subfactor,output)
            tracts = load(tract_file);
            indices = 1:subfactor:length(tracts.FList);
            tracts.TractAng = tracts.TractAng(indices);
            tracts.TractFA = tracts.TractFA(indices);
            tracts.TractFE = tracts.TractFE(indices);
            tracts.TractGEO = tracts.TractGEO(indices);
            tracts.TractL = tracts.TractL(indices);
            tracts.TractLambdas = tracts.TractLambdas(indices);
            tracts.TractMD = tracts.TractMD(indices);
            tracts.Tracts = tracts.Tracts(indices);
            try
                tracts.TractsFOD = tracts.TractsFOD(indices);
            catch
            end
            tracts.FList = 1:length(indices);
            save(output,'-struct','tracts');
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
function [bmat,LRKernel,HRKernel] = mDRLMT_MakeDKIKernel_multicomp(data,nreconstruction_vertices,lambdas,K,offset,isoDs,shell_data)

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
% lr_scheme = gen_scheme(min(length(data.bvals),90),4);
lr_scheme = gen_scheme(min(length(data.bvals),45),4);
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
            anglesFi = [phi(i), theta(i)];%*(180/pi); % in degrees
            Kernel{ij}(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K,offset);
        end
        for l=1:length(isoDs)
            Kernel{ij}(:,end+1) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
        end

        Kernel_LR{ij} = zeros(ndirections(ij),length(phi_LR));
        for i=1:length(phi_LR)
            anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
            Kernel_LR{ij}(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K,offset);
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
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K,offset);
    end
    for l=1:length(isoDs)
        HRKernel(:,length(phi)+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end

    for i=1:length(phi_LR)
        anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
        LRKernel(:,i) = create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K,offset);
    end
    for l=1:length(isoDs)
        LRKernel(:,length(phi_LR)+l) = create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
            bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
    end
end

end

% Private function assembling the NODDI based kernels.
function [bmat,LRKernel,HRKernel,super_scheme] = mDRLMT_MakeNODDIKernel_multicomp(data,nreconstruction_vertices,noddi_values,isoDs,shell_data)
global MRIToolkit
if(isempty(which('SynthMeasWatsonSHStickTortIsoV_B0')))
    if(isfield(MRIToolkit,'noddi_path') && exist(MRIToolkit.noddi_path,'dir') > 0)
        addpath(genpath(MRIToolkit.noddi_path))
    else
        error('Cannot find the NODDI toolbox. Please, add it to the MATLAB path');
    end
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
% Project:   High Angular Resolution Diffusion Imaging Tools
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
function [S, D] = create_signal_multi_tensor_dki (ang, f, Eigenvalues, b, grad, S0, SNR, add_noise, K, offset)
% -Normalizing the gradient vector and then transforming the b-value.
% -This part is only for non-unitary gradient vectors
% Transform_b_value_and_grad = repmat(sqrt(diag(grad*grad')+eps), [1 3]);
% grad = grad./Transform_b_value_and_grad;
% b = b.*(Transform_b_value_and_grad(:,1)).^2;
% grad = grad(:,[2 1 3]);
% S = 0;
% Nfibers = length(f);
% for i = 1:Nfibers
%     S = S + f(i)*EDTI_Library.SimulateSignalDTIIsotropicK(b,grad,Eigenvalues,K,flip([pi/2 0]+ang(i,:)));
% end
% D = 0;
% return
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
    S = S + f(i)*(exp(-b.*diag(grad*D*grad')+b2D2K)+offset);
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

if(N == 300)
    %     disp('Temp modification of gen_scheme');
    scheme.vert = load('dir300.txt');
    return
end

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

mult = sign(scheme.vert(:,3));
scheme.vert = scheme.vert .* mult;

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

% Taken from "hardi tools" of Erick Canales-Rodriguez
function fODF = ADT_deconv_RLdamp_1D_noEP_init(fODF, Signal, Kernel,KernelT,KTK, Niter,nu)
% fODF: purely anisotropic part

% KernelT = Kernel'; % only one time

fzero = single(1e-06);
KernelS = KernelT*Signal;
V = single(8);
% KTK = KernelT*Kernel;

mu = max(0, 1 - 4*std(Signal) );
nuV = nu^V;
last_odf = fODF;
my_eps = 1e-4*max(Signal);
for i = 1:Niter
    % --- approach: Flavio Dell’Acqua
    Dem = KTK*fODF;
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

% Taken from "hardi tools" of Erick Canales-Rodriguez
function fODF = ADT_deconv_RLdamp_1D_noEP_white(Signal, Kernel, Niter)
% fODF: purely anisotropic part

fODF0 = ones(size(Kernel,2),1);

fODF = fODF0/sum(fODF0);
KernelT = Kernel'; % only one time

fzero = 1e-06;
a = 0.9;
T = 0.15; N = 1/2;
for i = 1:Niter
    Reblurred = Kernel*fODF;
    % --- approach: L. White
    Rkk = ( -2/T )*(Signal.*log(Reblurred./(Signal + eps)) - Reblurred + Signal);
    R = min(1, Rkk);
    Uk = (R.^(N-1)).*(N - (N-1)*R);
    % ------
    K = 1 + Uk.*( (Signal - Reblurred)./(Reblurred + eps) );
    RL_factor = (KernelT*K).^a;
    fODF = max(fzero, fODF.*RL_factor); % positivity
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

% Using single precision caused some strange index shifting
% conn = repmat(single(0),[1 num_tracts + Total_P]);
conn = repmat(double(0),[1 num_tracts + Total_P]);

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

% A simple VTK read function to read VTK 4.2 and 5.1 pointsets produced by WMA
function [points_list,lines] = VTKImport_test(file_in)

fid = fopen(file_in,'rb');
str = fgets(fid);
parts = strsplit(strrep(str,newline,''),' ');
version = str2double(parts{end});
while(feof(fid) == false && contains(str,'POINTS') == false)
    str = fgets(fid);
end
nvert = sscanf(str,'%*s %d %*s', 1);
sp = ftell(fid);
% read vertices
A = single(fread(fid,nvert*3,'float','b'));
mp = ftell(fid);
while(feof(fid) == false && contains(str,'LINES') == false)
    str = fgets(fid);
end
% read lines
nlines = sscanf(str,'%*s %d %d');
if(version > 5)
    str = fgets(fid);
    LINES_INFO = fread(fid,nlines(1),'*int64','b');
    LINES_INFO = LINES_INFO;
    while(feof(fid) == false && contains(str,'CONNECTIVITY') == false)
        str = fgets(fid); % read datatype
    end
    lines = cell(nlines(1)-1,1);
    for line_id=1:length(lines)
        lines(line_id) = {(1+(LINES_INFO(line_id)):LINES_INFO(line_id+1)-1)}; % UPDATED
    end

    % Unused for now
    if(1 == 0)
        CONNECTIVITY_INFO = fread(fid,nlines(2),'*int64','b');
        while(feof(fid) == false && contains(str,'CELL_DATA') == false)
            str = fgets(fid); % read datatype
        end
        ncell = sscanf(str,'%*s %d');
        str = fgets(fid); % read FIELD
        str = fgets(fid); % read ClusterNumber
        nclusters = sscanf(str,'%*s %d %d %*s');
        cluster_type = char(sscanf(str,'%*s %*d %*d %s')');
        CLUSTER_INFO = fread(fid,nclusters(2),['*' cluster_type],'b');
        str = fgets(fid); % read empty line
        str = fgets(fid); % read EmbeddingColor
        ncolors = sscanf(str,'%*s %d %d %*s');
        COLORS_INFO = fread(fid,ncolors(1)*ncolors(2),'*uint8','b');
        str = fgets(fid); % read empty line
        str = fgets(fid); % read EmbeddingCoordinate
        ncoordinates = sscanf(str,'%*s %d %d %*s');
        COORD_INFO = fread(fid,ncoordinates(1)*ncoordinates(2),'*float','b');
        CELL_INFO = fread(fid,ncell,'*int64','b');
    end
else
    lines = cell(1,1);
    cnt = 1;
    while(cnt <= nlines)
        K = fread(fid,1,'int','b');
        if(K == 0)
            continue
        end
        lines(cnt) = {fread(fid,K,'*int','b')};
        cnt = cnt+1;
    end
end

fclose(fid);

A = reshape(A,3,nvert);
points_list = A';
end

% Ensure a vector is in column format
function col = ascolumn(in_vector)
if(~ismatrix(in_vector))
    warning('ascolumn works only with vectors - no nd matrices');
    col = [];
    return
end
if(size(in_vector,2) > size(in_vector,1))
    col = in_vector';
else
    col = in_vector;
end
end

% Tensor 2 Symmetric Tensor
function Dtensor = D2Dtensor(D)
Dtensor = [D(1) D(2)/2 D(3)/2
    D(2)/2 D(4) D(5)/2
    D(3)/2 D(5)/2 D(6)];
end

% Utility function to derive DTI metrics
function [MD,FA,DEC,AD,RD,LAMBDAS,L1] = DW_ComputeTensorMetrics(p,tensor_indexes)
sizes = size(p);
ndims = length(sizes);

rows = 1;
for j=1:ndims-1
    rows = rows*sizes(j);
end

newp = reshape(p,rows,sizes(end));
newp = newp(:,tensor_indexes);

MD = zeros(rows,1);
FA = zeros(rows,1);
if(nargout > 2)
    DEC = zeros(rows,3);
end

if(nargout > 3)
    AD = zeros(rows,1);
    RD = zeros(rows,1);
    LAMBDAS = zeros(rows,3);
    L1 = zeros(rows,3);
end

for j=1:size(newp,1)
    D = [newp(j,1) newp(j,2)/2 newp(j,3)/2
        newp(j,2)/2 newp(j,4) newp(j,5)/2
        newp(j,3)/2 newp(j,5)/2 newp(j,6)];
    %         D = [newp(j,1) newp(j,2) newp(j,3)
    %         newp(j,2) newp(j,4) newp(j,5)
    %         newp(j,3) newp(j,5) newp(j,6)];
    if(sum(~isfinite(D(:))) > 0 || sum(D(:)) == 0)
        continue
    end
    [autovett,autoval] = eig(D);
    autoval = diag(autoval);
    if(sum(autoval<0) == length(autoval))
        autoval = -autoval;
        autovett = -autovett;
    end
    autoval(autoval<0)=0;
    MD(j) = mean(autoval);
    FA(j) = sqrt(1.5*((autoval(1)-MD(j))^2+(autoval(2)-MD(j))^2+(autoval(3)-MD(j))^2)/(sum(autoval.^2)));
    %         FA(j) = sqrt(0.5*((autoval(1)-autoval(2))^2+(autoval(2)-autoval(3))^2+(autoval(1)-autoval(3))^2)/(sum(autoval.^2)));
    if(nargout > 2)
        [~,IX] = sort(autoval,'descend');
        DEC(j,:) = autovett(:,IX(1))*FA(j);
        if(nargout > 3)
            AD(j) = autoval(IX(1));
            RD(j) = 0.5*sum(autoval(IX(2:3)));
            LAMBDAS(j,:) = autoval(IX);
            L1(j,:) = autovett(IX(1),:);
        end
    end
end

MD = reshape(MD,sizes(1:end-1));
FA = reshape(FA,sizes(1:end-1));
if(nargout > 2)
    DEC = reshape(DEC,[sizes(1:end-1) 3]);
    if(nargout > 3)
        AD = reshape(AD,sizes(1:end-1));
        RD = reshape(RD,sizes(1:end-1));
        LAMBDAS = reshape(LAMBDAS,[sizes(1:end-1) 3]);
        L1 = reshape(L1,[sizes(1:end-1) 3]);
    end
end
end

% An L2 regularized NNLS function
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

% Compute a density weighted map from a set of streamlines
function E_DTI_DensityMaps(base_file,tract_file,out_file)
% filt_type 0: pass once in params.mask1
% filt_type 1: pass at least twice in params.mask1
% filt_type 2: start in params.mask1 and end in params.mask2

load(tract_file,'Tracts');
load(base_file,'FA','VDims');

[sx,sy,sz] = size(FA);
DensityMap = zeros([sx sy sz]);

try

    for tid=1:length(Tracts)
        tract = Tracts{tid};
        for l=1:size(tract,1)
            point = round(tract(l,:)./VDims);
            DensityMap(point(1),point(2),point(3)) = DensityMap(point(1),point(2),point(3)) + 1;
        end
    end

    DensityMap = DensityMap/max(DensityMap(:));
    DM.img = DensityMap;
    DM.VD = VDims;
    EDTI.WriteNifti(DM,out_file);

catch err
    disp(err.message);
end

end

function out = my_help(fname)
if(isdeployed)
    out = fname;
else
    out = help(fname);
end
end

% From Sparse-Wars
function fODF = sph_deconv_tv_motor_all(Signal, Kernel, fODF0, Niter, ncoils, MRI_recon_type)

% Erick J Canales-Rodríguez, Alessandro Daducci, Stamatios N Sotiropoulos, Emmanuel Caruyer, Santiago Aja-Fernández, Joaquim Radua, Yasser Iturria-Medina, Lester Melie-García, Yasser Alemán-Gómez, Jean-Philippe Thiran, Salvador Sarró, Edith Pomarol-Clotet, Raymond Salvador. Spherical Deconvolution of Multichannel Diffusion MRI Data with Non-Gaussian Noise Models and Spatial Regularization. PLoS ONE, 2015, 10(10): e0138910.

% ----------------------------------------------------
% fODF_4D: purely anisotropic part

if ( strcmp( MRI_recon_type, 'SMF-SENSE-based') ) || (ncoils == 1)
    n_order = 1;
elseif strcmp( MRI_recon_type, 'SoS-GRAPPA-based')
    n_order = ncoils;
end

fODF = fODF0;
Reblurred = Kernel*fODF;

fzero = 0;
KernelT = Kernel'; % only one time

lambda_aux = 0;

%%________________________ Main Algorithm ______________________________ %%
sigma0 = 1/15;
sigma2 = sigma0^2;
lambda = sigma2;
epsilon = 1e-7;
Data_2d = Signal;
N = size(Data_2d, 1);
sigma2 = sigma2*ones(size(Data_2d));
% --------------
Reblurred_S = (Signal.*Reblurred)./sigma2;
for i = 1:Niter
    %     display(['Iter -> ' num2str(i) ' of ' num2str(Niter)]);
    % ------------------- R-L deconvolution part -------------------- %

    fODFi = fODF;

    Ratio = mBessel_ratio(n_order,Reblurred_S);

    RL_factor = KernelT*( Data_2d.*( Ratio ) )./( KernelT*(Reblurred) + eps);

    fODF = fODFi.*RL_factor;

    fODF = max(fzero, fODF); % positivity

    %if i <= 100
    %cte = sqrt(1 + 2*n_order*mean(sigma2(:)));
    %cte = sqrt(1 + n_order*mean(sigma2(:)));
    %      cte = 1;
    %      fODF =  cte*fODF./repmat( sum(fODF,1) + eps, [size(fODF,1), 1] );
    % Energy preservation at each step, which included the bias on the S0 image.
    % This step is used only during the first iterations to stabilize
    % the recovery at the begining...
    % end

    % --------------------- Update of variables --------------------- %
    Reblurred = Kernel*fODF;
    Reblurred_S = (Data_2d.*Reblurred)./sigma2;

    % ---------------- noise ---------------- %
    sigma2_i = (1/N)*sum( (Data_2d.^2 + Reblurred.^2)/2 - (sigma2.*Reblurred_S).*Ratio, 1)./n_order;
    sigma2_i = min((1/8)^2, max(sigma2_i,(1/80)^2)); % estimator in the interval sigma = [1/SNR_min, 1/SNR_max],
    % where SNR_min = 8 and SNR_max = 80
    % --------------------------------------- %
    sigma2 = repmat(sigma2_i,[size(Signal,1), 1]);

end

fODF = fODF./repmat( sum(fODF,1) + eps, [size(fODF,1), 1] ); % energy preservation

%%____________________ End of Main Algorithm ___________________________ %%

end

% From Sparse Wars
function [fxyz] = grad(M)
% grad - gradient operator
fx = M([2:end end],:,:)-M;
fy = M(:,[2:end end],:)-M;
fz = M(:,:,[2:end end])-M;

fxyz = cat(4,fx,fy,fz);
end

% From Sparse Wars
function fd = div(P)
% div - divergence operator

Px = squeeze(P(:,:,:,1));
Py = squeeze(P(:,:,:,2));
Pz = squeeze(P(:,:,:,3));

fx = Px-Px([1 1:end-1],:,:);
fx(1,:,:)   = Px(1,:,:);    % boundary
fx(end,:,:) = -Px(end-1,:,:);
% ---
fy = Py-Py(:,[1 1:end-1],:);
fy(:,1,:)   = Py(:,1,:);    % boundary
fy(:,end,:) = -Py(:,end-1,:);
% ---
fz = Pz-Pz(:,:,[1 1:end-1]);
fz(:,:,1)   = Pz(:,:,1);    % boundary
fz(:,:,end) = -Pz(:,:,end-1);
% ---
fd = fx+fy+fz;
end

% from Sparse Wars
function y = mBessel_ratio(n,x)
% y = mBessel{n}(x)/mBessel{n-1}(x) = besseli(n,x)./besseli(n-1,x)
% Fast evaluation using the Perron's Continued Fraction equation.
% For more details see: http://www.ams.org/journals/mcom/1978-32-143/S0025-5718-1978-0470267-9/

y = x./( (2*n + x) - ( 2*x.*(n+1/2)./ ( 2*n + 1 + 2*x - ( 2*x.*(n+3/2)./ ( 2*n + 2 + 2*x - ( 2*x.*(n+5/2)./ ( 2*n + 3 + 2*x ) ) ) ) ) ) );
end

