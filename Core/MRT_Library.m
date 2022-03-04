%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

% This class transforms the processing methods originally implemented in ExploreDTI
% into a library, made consistent into a class enforcing its consistency,
% ease of use and availability for scripting / command line tools without
% need for a GUI.
% Author: Alberto De Luca - alberto@isi.uu.nl - alberto.deluca.06@gmail.com
% Many of the methods here included were originally implemented by Alexander Leemans.
% First version 29/10/2019 as part of the class EDTI
% 24/10/2020: All ExploreDTI are adapted and moved to this class as static
% methods
classdef MRT_Library < handle
    
    methods(Static)
        
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
            
            super_scheme = MRT_Library.gen_scheme(nreconstruction_vertices,4); % the reconstruction scheme. Change 300 to any number
            HRKernel = zeros(sum(ndirections),nreconstruction_vertices+length(isoDs));
            [phi, theta] = cart2sph(super_scheme.vert(:,1),super_scheme.vert(:,2),super_scheme.vert(:,3)); % polar decomposition
            lr_scheme = MRT_Library.gen_scheme(min(length(data.bvals),90),4);
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
                        Kernel{ij}(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                            bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K);
                    end
                    for l=1:length(isoDs)
                        Kernel{ij}(:,end+1) = MRT_Library.create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                            bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
                    end
                    
                    Kernel_LR{ij} = zeros(ndirections(ij),length(phi_LR));
                    for i=1:length(phi_LR)
                        anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
                        Kernel_LR{ij}(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                            bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K);
                    end
                    for l=1:length(isoDs)
                        Kernel_LR{ij}(:,end+1) = MRT_Library.create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
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
                    HRKernel(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                        bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K);
                end
                for l=1:length(isoDs)
                    HRKernel(:,length(phi)+l) = MRT_Library.create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                        bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
                end
                
                for i=1:length(phi_LR)
                    anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
                    LRKernel(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                        bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K);
                end
                for l=1:length(isoDs)
                    LRKernel(:,length(phi_LR)+l) = MRT_Library.create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
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
                R = MRT_Library.RotMatrix(phi(i),theta(i));
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
                R = MRT_Library.RotMatrix(phi(i),theta(i));
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
        
        function noisy_signal = AddRicianNoise(S,sigma)
            
            standar_deviation = sigma.*(ones(length(S),1));
            med = zeros(length(S),1);
            
            er1 = normrnd(med, standar_deviation);
            er2 = normrnd(med, standar_deviation);
            noisy_signal = sqrt((S + er1).^2 + er2.^2);
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
            
            if(N == 300)
%                 disp('Temp modification of gen_scheme');
                scheme.vert = load('dir300.txt');
                return
            end
            
            if size(N,1) == 1 & size(N,2) == 1
                P = c2s(MRT_Library.equidistribute(N));
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
                scheme.sh = [ scheme.sh MRT_Library.eval_SH(l, scheme.el, scheme.az)' ];
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
            
            s = MRT_Library.eval_ALP(l, el).*s';
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
        
        % Calculate the dampening factor for Richardson Lucy, inspired by
        % discussion with Flavio Dell'Acqua (thanks!)
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
        
        function dec = ConstrainedDeconvolution(X,S,lambda)
            iters = 30;
            cost = NaN(iters,1);

            Xp = X;
            Sp = S;

            REG = 0.001*lambda*ones(size(Xp,2),1);
            XT = Xp';
            XX = XT*Xp;
            XS = XT*Sp;
            for iter=1:iters   
               p = inv(XX + diag(REG))*XS;   
               cost(iter) = sum((S-Xp*p).^2);
%                disp(cost(iter))

               neg = p < eps;
               REG(neg) = lambda;
               REG(~neg) = 0.001*lambda;
               
               if(iter > 10 && cost(iter) == cost(iter-1))
                   break
               end
            end
            dec = p;
            dec(dec<eps) = 0;
        end
        
        function out = my_help(fname)
            if(isdeployed)
                out = fname;
            else
                out = help(fname);
            end
        end
        
    end
    
end
