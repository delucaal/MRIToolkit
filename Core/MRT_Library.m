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
            
            super_scheme = MRT_Library.gen_scheme(nreconstruction_vertices,4); % the reconstruction scheme. Change 300 to any number
            HRKernel = zeros(sum(ndirections),nreconstruction_vertices+length(isoDs));
            [phi, theta] = cart2sph(super_scheme.vert(:,1),super_scheme.vert(:,2),super_scheme.vert(:,3)); % polar decomposition
            % lr_scheme = gen_scheme(min(length(data.bvals),90),4);
            lr_scheme = MRT_Library.gen_scheme(min(length(data.bvals),45),4);
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
                        Kernel{ij}(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                            bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K,offset);
                    end
                    for l=1:length(isoDs)
                        Kernel{ij}(:,end+1) = MRT_Library.create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                            bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0);
                    end
                    
                    Kernel_LR{ij} = zeros(ndirections(ij),length(phi_LR));
                    for i=1:length(phi_LR)
                        anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
                        Kernel_LR{ij}(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                            bmat{ij}(:,4), bmat{ij}(:,1:3), 1, 0, 0, K,offset);
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
                        bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K,offset);
                end
                for l=1:length(isoDs)
                    HRKernel(:,length(phi)+l) = MRT_Library.create_signal_multi_tensor([0 0], 1, [isoDs(l) isoDs(l) isoDs(l)], ...
                        bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0);
                end
                
                for i=1:length(phi_LR)
                    anglesFi = [phi_LR(i), theta_LR(i)]*(180/pi); % in degrees
                    LRKernel(:,i) = MRT_Library.create_signal_multi_tensor_dki(anglesFi, 1, lambdas, ...
                        bmat{1}(:,4), bmat{1}(:,1:3), 1, 0, 0, K,offset);
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
            if(exist('offset','var') < 1)
                offset = 0;
            end

            S = 0;
            Nfibers = length(f);
            f = f/sum(f);
            b2D2K = 1/6*b.^2.*mean(Eigenvalues.^2)*K;
            for i = 1:Nfibers
                phi(i) = ang(i, 1);
                theta(i) = ang(i, 2);
                R = MRT_Library.RotMatrix(phi(i),theta(i));
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
            my_eps = 1e-4*max(fODF);%max(Signal);
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

        % Remove known extensions from a filename to only retain its prefix
        function [prefix,extension] = FilenameNoExtension(filename)
            known_extensions = {'.nii','.gz','.vtk','.trk'};
            extension = [];
            prefix = filename;
            for ext=1:length(known_extensions)
                if(contains(filename,known_extensions{ext}))
                    extensions = cat(2,extensions,known_extensions{ext});
                    prefix = strrep(prefix,known_extensions{ext},'');
                end
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
        
        % An own version of constrained (non-negative) least squares 
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
        
        % Perform fiber tractography of 2 SH-based FODs (mFOD)
        function WholeBrain_mFOD_Tractography(reference_mat,fod_1,fod_2,p,f_out)
            tic
            
            f_in = reference_mat;
            
            load(f_in,'VDims')
            
            if(ischar(fod_1))
                CSD_FOD = EDTI_Library.E_DTI_read_nifti_file(fod_1);
            else
                CSD_FOD = fod_1;
                clear fod_1;
            end

            if(ischar(fod_2))
                CSD_FOD2 = EDTI_Library.E_DTI_read_nifti_file(fod_2);
            else
                CSD_FOD2 = fod_2;
                clear fod_2;
            end
            
            SR = p.SeedPointRes;
            mask_s = ~isnan(CSD_FOD(:,:,:,1)) | ~isnan(CSD_FOD(:,:,:,2));
            
            disp('Calculating seed points...')
            if(isfield(p,'SeedMask') && ~isempty(p.SeedMask))
                if(ischar(p.SeedMask))
                    seed_mask = EDTI_Library.E_DTI_read_nifti_file(p.SeedMask);
                    p.SeedMask = seed_mask;
                end
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(p.SeedMask > 0, SR, VDims, p);
            else
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(mask_s, SR, VDims, p);
            end
            
            if isempty(seedPoint)
                disp('No seed points found - processing stopped.')
            else
                disp('Seed point calculation done.')
            end
            
            disp('Calculating trajectories...')
            
            stepSize = p.StepSize;
            threshold = p.blob_T;
            maxAngle = p.AngleThresh;
            lengthRange = [p.FiberLengthRange(1) p.FiberLengthRange(2)];
            v2w = diag(VDims); v2w(4,4) = 1;
            
            t = MultiFieldSHTracker(v2w);

            t.setData(CSD_FOD,CSD_FOD2);
            t.setParameters(stepSize, threshold, maxAngle, lengthRange); t.setProgressbar(false);
            [Tracts, TractsFOD] = t.track(seedPoint);
            
            num_tracts = size(Tracts,2);
            TractL = cell(1,num_tracts);
            for i = 1:num_tracts
                TractL{i} = single(size(Tracts{i},1));
            end
            
            TL = cell2mat(TractL(:))*stepSize;
            IND = and(TL>=lengthRange(1),TL<=lengthRange(2));
            Tracts = Tracts(IND);
            TractsFOD = TractsFOD(IND);
            
            num_tracts = size(Tracts,2);
            TractL = cell(1,num_tracts);
            for i = 1:num_tracts
                TractL{i} = single(size(Tracts{i},1));
            end
            
            disp('Trajectory calculations done.')
            
            disp('Processing diffusion info along the tracts...');
            
            load(f_in,'DT','MDims');
            
            [TractFA, TractFE, TractAng, TractGEO, TractLambdas, TractMD] =...
                EDTI_Library.E_DTI_diff_measures_vectorized(Tracts, VDims, DT);
            
            % TractL = cellfun(@length, Tracts, 'UniformOutput', 0);
            
            FList = (1:length(Tracts))';
            
            TractMask = repmat(0,MDims);
            
            for i = 1:length(FList)
                IND = unique(sub2ind(MDims,...
                    round(double(Tracts{i}(:,1))./VDims(1)),...
                    round(double(Tracts{i}(:,2))./VDims(2)),...
                    round(double(Tracts{i}(:,3))./VDims(3))));
                TractMask(IND) = TractMask(IND) + 1;
            end
            
            disp('Processing diffusion info along the tracts done.');
            
            Parameters.Step_size = stepSize;
            Parameters.FOD_threshold = threshold;
            Parameters.Angle_threshold = maxAngle;
            Parameters.Length_range = lengthRange;
            Parameters.Seed_resolution = SR;
            
            
            disp('Saving trajectories...')
            save(f_out,'FList','Tracts','TractL','TractFE',...
                'TractFA','TractAng','TractGEO','TractMD',...
                'TractLambdas','TractMask','VDims','TractsFOD','Parameters','-v7.3');
            disp('Saving trajectories done.')
            
            ti = toc;
            
            m=ti/60;
            disp(['Tracking (CSD - FOD interpolation) computation time was ' num2str(m) ' min.'])
            
        end

        % Originally from ExploreDTI: Perform fiber tractography of any FOD in SH - Adapted
        % and parallelized. Supports multiple trackers
        function WholeBrainFODTractography_par(reference_mat,CSD_FOD_or_file,p,f_out)
            global MRIToolkit
            
            tic
            
            f_in = reference_mat;
            
            load(f_in,'VDims')
            
            if(ischar(CSD_FOD_or_file))
                CSD_FOD = EDTI_Library.E_DTI_read_nifti_file(CSD_FOD_or_file);
            else
                CSD_FOD = CSD_FOD_or_file;
                clear CSD_FOD_or_file;
            end
            
            SR = p.SeedPointRes;
            mask_s = ~isnan(CSD_FOD(:,:,:,1));
            
            disp('Calculating seed points...')
            if(isfield(p,'SeedMask') && ~isempty(p.SeedMask))
                if(ischar(p.SeedMask))
                    seed_mask = EDTI_Library.E_DTI_read_nifti_file(p.SeedMask);
                    p.SeedMask = seed_mask;
                end
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(p.SeedMask > 0, SR, VDims, p);
            else
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(mask_s, SR, VDims, p);
            end
            
            if isempty(seedPoint)
                disp('No seed points found - processing stopped.')
            else
                disp('Seed point calculation done.')
            end
            
            disp('Calculating trajectories...')
            
            stepSize = p.StepSize;
            threshold = p.blob_T;
            maxAngle = p.AngleThresh;
            lengthRange = [p.FiberLengthRange(1) p.FiberLengthRange(2)];
            v2w = diag(VDims); v2w(4,4) = 1;
            
            pool = gcp;
            if(isempty(pool))
                % Could not create parallel pool, resorting to original
                % implementation
                EDTI_Library.WholeBrainFODTractography_par(reference_mat,CSD_FOD_or_file,p,f_out);
                return
            else
                % Divide seeds accordingly
                Step = ceil(size(seedPoint,2)/pool.NumWorkers);
                seedPointsSplit = cell(pool.NumWorkers,1);
                for ij=1:pool.NumWorkers
                    seedPointsSplit(ij) = {seedPoint(:,(ij-1)*Step+1:min(ij*Step,size(seedPoint,2)))};
                end

                Tracts = cell(pool.NumWorkers,1);
                TractsFOD = cell(pool.NumWorkers,1);
                parfor block = 1:pool.NumWorkers
                    t = [];
                    if(isfield(MRIToolkit,'fibertracker') < 1 || isfield(MRIToolkit.fibertracker,'type') < 1 ...
                            || isempty(MRIToolkit.fibertracker.type))
                        t = SHTracker(v2w);
                    else
%                         disp('Not yet supported')
%                         continue
                          disp('You are using a very experimental fiber tracker!');
                          t = DistProbSHTracker(v2w);
                          t.setNumberOfIterations(10);
                          t.setParametersSD(0.3,10);
                          t.weight_mode = 0;

%                         eval(['t = ' MRIToolkit.fibertracker.type '(v2w);']);
%                         if(isfield(MRIToolkit.fibertracker,'parameters'))
%                             fields = fieldnames(MRIToolkit.fibertracker.parameters);
%                             for field_id=1:length(fields)
%                                 eval(['t.' fields{field_id} '= MRIToolkit.fibertracker.parameters.' fields{field_id} ';']);
%                             end
%                         end
                    end
                    t.setData(CSD_FOD);
                    t.setParameters(stepSize, threshold, maxAngle, lengthRange); t.setProgressbar(false);
                    [LTracts, LTractsFOD] = t.track(seedPointsSplit{block});
                    Tracts(block) = {LTracts};
                    TractsFOD(block) = {LTractsFOD};
                end

                for ix=2:length(Tracts)
                    Tracts(1) = {cat(2,Tracts{1},Tracts{ix})};
                    TractsFOD(1) = {cat(2,TractsFOD{1},TractsFOD{ix})};
                end

                Tracts = Tracts{1};
                TractsFOD = TractsFOD{1};

                num_tracts = size(Tracts,2);
                TractL = cell(1,num_tracts);
                for i = 1:num_tracts
                    TractL{i} = single(size(Tracts{i},1));
                end
                
                TL = cell2mat(TractL(:))*stepSize;
                IND = and(TL>=lengthRange(1),TL<=lengthRange(2));
                Tracts = Tracts(IND);
                TractsFOD = TractsFOD(IND);
                
                num_tracts = size(Tracts,2);
                TractL = cell(1,num_tracts);
                for i = 1:num_tracts
                    TractL{i} = single(size(Tracts{i},1));
                end
            end

            disp('Trajectory calculations done.')
            
            disp('Processing diffusion info along the tracts...');
            
            load(f_in,'DT','MDims');
            
            [TractFA, TractFE, TractAng, TractGEO, TractLambdas, TractMD] =...
                EDTI_Library.E_DTI_diff_measures_vectorized(Tracts, VDims, DT);
            
            % TractL = cellfun(@length, Tracts, 'UniformOutput', 0);
            
            FList = (1:length(Tracts))';
            
            TractMask = repmat(0,MDims);
            
            for i = 1:length(FList)
                IND = unique(sub2ind(MDims,...
                    round(double(Tracts{i}(:,1))./VDims(1)),...
                    round(double(Tracts{i}(:,2))./VDims(2)),...
                    round(double(Tracts{i}(:,3))./VDims(3))));
                TractMask(IND) = TractMask(IND) + 1;
            end
            
            disp('Processing diffusion info along the tracts done.');
            
            Parameters.Step_size = stepSize;
            Parameters.FOD_threshold = threshold;
            Parameters.Angle_threshold = maxAngle;
            Parameters.Length_range = lengthRange;
            Parameters.Seed_resolution = SR;
            
            
            disp('Saving trajectories...')
            save(f_out,'FList','Tracts','TractL','TractFE',...
                'TractFA','TractAng','TractGEO','TractMD',...
                'TractLambdas','TractMask','VDims','TractsFOD','Parameters','-v7.3');
            disp('Saving trajectories done.')
            
            ti = toc;
            
            m=ti/60;
            disp(['Tracking (CSD - FOD interpolation) computation time was ' num2str(m) ' min.'])
            
        end
        
        % Originally from ExploreDTI: Perform fiber tractography of any FOD in SH - Adapted
        % and parallelized. Supports multiple trackers
        function WholeBrainFODProbTractography_par(reference_mat,CSD_FOD_or_file,p,f_out)
            global MRIToolkit
            
            tic
            
            f_in = reference_mat;
            
            load(f_in,'VDims')
            
            if(ischar(CSD_FOD_or_file))
                CSD_FOD = EDTI_Library.E_DTI_read_nifti_file(CSD_FOD_or_file);
            else
                CSD_FOD = CSD_FOD_or_file;
                clear CSD_FOD_or_file;
            end
            
            SR = p.SeedPointRes;
            mask_s = ~isnan(CSD_FOD(:,:,:,1));
            
            disp('Calculating seed points...')
            if(isfield(p,'SeedMask') && ~isempty(p.SeedMask))
                if(ischar(p.SeedMask))
                    seed_mask = EDTI_Library.E_DTI_read_nifti_file(p.SeedMask);
                    p.SeedMask = seed_mask;
                end
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(p.SeedMask > 0, SR, VDims, p);
            else
                seedPoint = EDTI_Library.E_DTI_Get_Seeds_WBT(mask_s, SR, VDims, p);
            end
            
            if isempty(seedPoint)
                disp('No seed points found - processing stopped.')
            else
                disp('Seed point calculation done.')
            end
            
            disp('Calculating trajectories...')
            
            iterations = 10;
            stepSize = p.StepSize;
            threshold = p.blob_T;
            maxAngle = p.AngleThresh;
            lengthRange = [p.FiberLengthRange(1) p.FiberLengthRange(2)];
            v2w = diag(VDims); v2w(4,4) = 1;
            
            pool = gcp;
            if(isempty(pool))
                % Could not create parallel pool, resorting to original
                % implementation
                EDTI_Library.WholeBrainFODTractography_par(reference_mat,CSD_FOD_or_file,p,f_out);
                return
            else
                % Divide seeds accordingly
                Step = ceil(size(seedPoint,2)/pool.NumWorkers);
                seedPointsSplit = cell(pool.NumWorkers,1);
                for ij=1:pool.NumWorkers
                    seedPointsSplit(ij) = {seedPoint(:,(ij-1)*Step+1:min(ij*Step,size(seedPoint,2)))};
                end

                mat_files_to_join = {};
                for iters=1:iterations
                    disp(['Iteration ' num2str(iters)]);
                    Tracts = cell(pool.NumWorkers,1);
                    TractsFOD = cell(pool.NumWorkers,1);
                    parfor block = 1:pool.NumWorkers
                        t = [];
                        if(isfield(MRIToolkit,'fibertracker') < 1 || isfield(MRIToolkit.fibertracker,'type') < 1 ...
                                || isempty(MRIToolkit.fibertracker.type))
                            t = DistProbSHTracker(v2w);
    %                         t = MultiPeakTracker(v2w);
                        else
                            disp('Not yet supported')
                            continue
    %                         eval(['t = ' MRIToolkit.fibertracker.type '(v2w);']);
    %                         if(isfield(MRIToolkit.fibertracker,'parameters'))
    %                             fields = fieldnames(MRIToolkit.fibertracker.parameters);
    %                             for field_id=1:length(fields)
    %                                 eval(['t.' fields{field_id} '= MRIToolkit.fibertracker.parameters.' fields{field_id} ';']);
    %                             end
    %                         end
                        end
                        t.setData(CSD_FOD);
                        t.setParameters(stepSize, threshold, maxAngle, lengthRange); t.setProgressbar(false);
                        t.setNumberOfIterations(1);
    %                     t.setParametersSD(stepSize/10,maxAngle/5);
                        t.setParametersSD(stepSize/5,maxAngle/10);
                        [LTracts, LTractsFOD] = t.track(seedPointsSplit{block});
                        Tracts(block) = {LTracts};
                        TractsFOD(block) = {LTractsFOD};
                    end
    
                    for ix=2:length(Tracts)
                        Tracts(1) = {cat(2,Tracts{1},Tracts{ix})};
                        TractsFOD(1) = {cat(2,TractsFOD{1},TractsFOD{ix})};
                    end
    
                    Tracts = Tracts{1};
                    TractsFOD = TractsFOD{1};
    
                    num_tracts = size(Tracts,2);
                    TractL = cell(1,num_tracts);
                    for i = 1:num_tracts
                        TractL{i} = single(size(Tracts{i},1));
                    end
                    
                    TL = cell2mat(TractL(:))*stepSize;
                    IND = and(TL>=lengthRange(1),TL<=lengthRange(2));
                    Tracts = Tracts(IND);
                    TractsFOD = TractsFOD(IND);
                    
                    num_tracts = size(Tracts,2);
                    TractL = cell(1,num_tracts);
                    for i = 1:num_tracts
                        TractL{i} = single(size(Tracts{i},1));
                    end
    
                    disp('Trajectory calculations done.')
                    
                    disp('Processing diffusion info along the tracts...');
                    
                    load(f_in,'DT','MDims');
                    
                    [TractFA, TractFE, TractAng, TractGEO, TractLambdas, TractMD] =...
                        EDTI_Library.E_DTI_diff_measures_vectorized(Tracts, VDims, DT);
                    
                    % TractL = cellfun(@length, Tracts, 'UniformOutput', 0);
                    
                    FList = (1:length(Tracts))';
                    
                    TractMask = repmat(0,MDims);
                    
                    for i = 1:length(FList)
                        IND = unique(sub2ind(MDims,...
                            round(double(Tracts{i}(:,1))./VDims(1)),...
                            round(double(Tracts{i}(:,2))./VDims(2)),...
                            round(double(Tracts{i}(:,3))./VDims(3))));
                        TractMask(IND) = TractMask(IND) + 1;
                    end
                    
                    disp('Processing diffusion info along the tracts done.');
                    
                    Parameters.Step_size = stepSize;
                    Parameters.FOD_threshold = threshold;
                    Parameters.Angle_threshold = maxAngle;
                    Parameters.Length_range = lengthRange;
                    Parameters.Seed_resolution = SR;
                    
                    
                    disp('Saving trajectories...')
                    mat_files_to_join = cat(2,mat_files_to_join,['mat_file' num2str(iters)]);
                    mat_files_to_join = cat(2,mat_files_to_join,strrep(f_out,'.mat',sprintf('_block%03d.mat',iters)));
                    save(strrep(f_out,'.mat',sprintf('_block%03d.mat',iters)),'FList','Tracts','TractL','TractFE',...
                        'TractFA','TractAng','TractGEO','TractMD',...
                        'TractLambdas','TractMask','VDims','TractsFOD','Parameters','-v7.3');
                    disp('Saving trajectories done.')
                    disp('End iteration');
                end
                ti = toc;
                
                m=ti/60;
                disp(['Tracking (CSD - FOD interpolation) computation time was ' num2str(m) ' min.'])
            
                disp('Now joining all tracts in a single file')
                mat_files_to_join = cat(2,mat_files_to_join,'output');
                mat_files_to_join = cat(2,mat_files_to_join,f_out);
                MRTTrack.ConcatenateMATTractographies(mat_files_to_join{:})
            end
        end        

        % Crop data given a specific bounding box
        function cropped_data = CropDataWithBoundingBox(Vol,extrapad)
           if(nargin < 2)
               extrapad = 0;
           end
            
           ref = Vol(:,:,:,1);
           thr = 0.01*prctile(ref(:),99);
           p1 = sum(sum(ref,2),3);
           p2 = squeeze(sum(sum(ref,1),3));
           p3 = squeeze(sum(sum(ref,1),2));
           
           mx = find(p1 > thr);
           my = find(p2 > thr);
           mz = find(p3 > thr);
           
           rx = max(1,mx(1)-extrapad):min(size(Vol,1),mx(end)+extrapad);
           ry = max(1,my(1)-extrapad):min(size(Vol,2),my(end)+extrapad);
           rz = max(1,mz(1)-extrapad):min(size(Vol,3),mz(end)+extrapad);
           
           cropped_data = Vol(rx,ry,rz,:);
           
        end

        % Adapted from ExploreDTI: helper function for Motion/Eddy/EPI correction
        % Performs the b-matrix rotation. Generalized to multiple rotations
        function [b, g, R] = Reorient_grad_and_b_matrix_rigid_rotation(b,g,Rot)
            R = eye(3);
            for x=1:length(Rot)
                M = zeros(4,3);
                M(1,:) = [Rot{x}{1} Rot{x}{2} Rot{x}{3}]*(180/pi);
                M(3,:) = 1;
                World_Trafo = EDTI_Library.E_DTI_Convert_aff_par_2_aff_transf_matrix(M);
                R = World_Trafo(1:3,1:3)*R;
            end

            b_old = b;
            g_old = g;
            
            for i=1:size(b,1)
                
                B = [b_old(i,1) b_old(i,2)/2 b_old(i,3)/2;...
                    b_old(i,2)/2 b_old(i,4) b_old(i,5)/2;...
                    b_old(i,3)/2 b_old(i,5)/2 b_old(i,6)];
                B_rot = R*B*(R');
                
                b(i,:) = [B_rot(1,1) 2*B_rot(1,2) 2*B_rot(1,3)...
                    B_rot(2,2) 2*B_rot(2,3) B_rot(3,3)];
                
            end
            
            
            for i=1:size(g,1)
                
                g(i,:) = g_old(i,:)*R';
                
            end
        end


        % From ExploreDTI: helper function for Motion/Eddy/EPI correction
        % This actually performs the EPI correction
        function data = TransformDataWithElastixParams(mat_file,out_name,ref_file,trafo_final_params)
            f_in = mat_file;
            
            [FOLD,FFO] = fileparts(f_in);
            
            par.temp_folder = fullfile(tempdir,['MRT_' num2str(randi(150)) '_' num2str(randi(150)) '_' num2str(randi(150))]);
            if(strcmp(par.temp_folder(end),'/'))
                par.temp_folder = par.temp_folder(1:end-1);
            end
            [FOL_FI,F] = fileparts(f_in);
            dir_temp = [par.temp_folder filesep 'Temp_' F];
            mkdir(dir_temp)

            par.dir_temp = dir_temp;
%             par.R2D.FN = ref_file;
            par.R2D.type = 1;
%             [suc, FN_nii] = EDTI_Library.E_DTI_SMECEPI_check_data_stuff(f_in,par);
            for_trafo = par;

%             if suc==0
%                 return;
%             end
            
%             fn_fixed = [for_trafo.dir_temp filesep 'Trafo_fixed.nii'];
%             fn_fixed_mask = [for_trafo.dir_temp filesep 'Trafo_fixed_mask.nii'];
%             [Fixed,VDims] = EDTI_Library.E_DTI_read_nifti_file(FN_nii);
%             VDims_E = VDims;
%             Fixed = single(Fixed);
%             Fixed = 10000*(Fixed/max(Fixed(:)));
%             EDTI_Library.E_DTI_write_nifti_file(Fixed,VDims,fn_fixed);
            
%             mask = single(Fixed>0);
%             EDTI_Library.E_DTI_write_nifti_file(mask,VDims,fn_fixed_mask);
            
            fn_moving = [for_trafo.dir_temp filesep 'Trafo_moving.nii'];
            
            LoadF = mat_file;
            
            [~,linked_files] = MRT_Library.ElastixParametersTree(trafo_final_params,1);

            Trafo_rig_result = {};
            for tpid=1:length(linked_files)
                params = MRT_Library.ReadElastixParameters(linked_files{tpid});
                ttype = MRT_Library.GetElastixParameter(params,'Transform');
                if(contains(ttype,'AffineDTITransform'))
                    Trafo_rig_result(end+1) = linked_files(tpid);
                end
            end
            
            Rig = cell(length(Trafo_rig_result));
            for ix=1:length(Rig)
                Q = textread(Trafo_rig_result{ix},'%s');
                Rig{ix}{1} = str2num(Q{7});
                Rig{ix}{2} = str2num(Q{6});
                Rig{ix}{3} = str2num(Q{8});
            end

            hdr = [];
            load(LoadF,'b','bval','info','g','NrB0','DWI','VDims','hdr','par','Mask_par')
        
%             if(~isempty(hdr))
%                 out.hdr = hdr;
%             end
            out=[];
            for ix=1:length(DWI)
                out.VD = VDims;
                out.img = DWI{ix};
                f2s = [dir_temp filesep 'Vol_' sprintf('%05d.nii',ix)];
                f2s_t = [dir_temp filesep 'Vol_' sprintf('%05d_trafo.nii',ix)];
                f2s_tp = [dir_temp filesep 'Vol_' sprintf('%05d_trafo_FP.nii',ix)]; 
                MRTQuant.WriteNifti(out,f2s);
                DW_Elastix_Transform(f2s,f2s_t,trafo_final_params);
                MRTQuant.ConformSpatialDimensions('nii_file',f2s_t,'output',f2s_tp)
                I = MRTQuant.LoadNifti(f2s_tp);
                VDims_E = I.VD;
                DWI(ix) = {I.img};
            end

            [b, g] = MRT_Library.Reorient_grad_and_b_matrix_rigid_rotation(b,g,Rig);
            
            if isnan(bval)
                diff_model=2;
            else
                diff_model=1;
            end
            
            [dummy, g] =  EDTI_Library.E_DTI_Get_Gradients_and_Bvalue(b, NrB0, diff_model);
            
            for i=1:length(DWI)
                if isempty(DWI{i})
                    disp('Errors encountered...')
                    EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                    suc = 0;
                    return;
                end
            end
            
            for i=1:length(DWI)
                DWI{i} = single(DWI{i});
                DWI{i}(DWI{i}==-1000)=nan;
            end
            
            
            dummy = false(size(DWI{1}));
            for i=1:length(DWI)
                dummy = or(dummy,isnan(DWI{i}));
                DWI{i}(isnan(DWI{i}))=0;
                DWI{i}(DWI{i}<0)=0;
            end

            par.cust_mask.NS = '';
            par.cust_mask.TS = '';
            par.mask_P = Mask_par;
            par.mask_P.NS.mfs = 5;
            par.mask_P.NS.NDWI = 0.7;
            par.mask_P.NS.DWI = 0.7;
            par.mask_P.TS.mfs = 5;
            par.mask_P.TS.NDWI = 0.7;
            par.mask_P.TS.DWI = 0.7;

            if isempty(par.cust_mask.TS)
                mask = EDTI_Library.E_DTI_Create_Mask_From_DWI_enhanced(DWI,NrB0,par.mask_P.TS.NDWI,par.mask_P.TS.DWI,par.mask_P.TS.mfs);
            else
                fn_cm = [FOLD filesep FFO par.cust_mask.TS];
                [mask, VDims, suc] = EDTI_Library.E_DTI_read_nifti_file(fn_cm);
                if suc==0
                    EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                    return;
                end
                mask = mask>0;
            end
            
            mask(dummy)=0;
            
            par_temp = par;
            %par_temp.TE = par.TE.TS;
            
            if diff_model==1
                [DT, DWIB0, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_from_DWI_b_mask(DWI,b,mask,par_temp,NrB0);
            elseif diff_model==2
                [DT, DWIB0, KT, outlier, chi_sq, chi_sq_iqr] = EDTI_Library.E_DTI_Get_DT_KT_from_DWI_b_mask_with_constraints(DWI,b,mask,g,par_temp,NrB0,VDims);
            end
            
            [FEFA, FA, FE, SE, eigval] = EDTI_Library.E_DTI_eigensystem_analytic(DT);
            
            g = g(NrB0+1:end,:);
            
            if ~isreal(FA)
                FA = real(FA);
            end
            if ~isreal(FEFA)
                FEFA = real(FEFA);
            end
            if ~isreal(FE)
                FE = real(FE);
            end
            if ~isreal(SE)
                SE = real(SE);
            end
            if ~isreal(eigval)
                eigval = real(eigval);
            end
            MDims = size(mask);
            VDims = VDims_E;
            
            f_out = out_name;
            par.Rig = Rig;
            
            max_DWI = 0;
            for i=1:length(DWI)
                max_DWI = max(max_DWI,max(DWI{i}(:)));
            end
            
            if max_DWI<=intmax('int16')
                for i=1:length(DWI)
                    DWI{i} = round(DWI{i});
                    DWI{i} = int16(DWI{i});
                end
            elseif max_DWI<=intmax('uint16')
                for i=1:length(DWI)
                    DWI{i} = round(DWI{i});
                    DWI{i}(DWI{i}<0)=0;
                    DWI{i} = uint16(DWI{i});
                end
            end
            
            
            try
                save(f_out,'DWI','VDims','b','bval','g','info','FEFA','NrB0','MDims',...
                    'FA','FE','SE','eigval','DT','outlier','DWIB0','chi_sq','chi_sq_iqr','par','hdr','Mask_par')
            catch me
                EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
                disp(me.message)
                return;
            end
            
            if diff_model==2
                save(f_out,'KT','-append')
            end
            
            EDTI_Library.E_DTI_remove_temp_f(for_trafo.dir_temp);
            if(nargout > 0)
                data = MRTQuant.EDTI_Data_2_MRIToolkit('mat_file',f_out,'do_preproc',0);
            end
        end

        % Concatenate multiple transformations of Elastix
        % Transformation_files should be in order: from first (initial step)
        % to last (final space)
        function start_index = ConcatenateElastixParams(transformation_files,output_folder,series_index,start_index)
            if(nargin < 3)
                start_index = 1;
                series_index = 1;
            end
            if(start_index == 1 && series_index == 1)
                % Create the output directory at the very beginning
                mkdir(output_folder);
            end
            if(iscell(transformation_files))
                % Recursively calls itself
                max_depth = zeros(length(transformation_files),1);
                for idx=1:length(transformation_files)
                    disp(['Recursion depth ' num2str(start_index) ' on ' transformation_files{idx}])
                    max_depth(idx) = MRT_Library.ConcatenateElastixParams(transformation_files{idx},output_folder,...
                        idx,1);
                end
                if(series_index == 1 && start_index == 1)
                    % Execute only once as the data has been copied
                    for idx=1:length(transformation_files)
                        for lev=1:max_depth(idx)
                            f2r = fullfile(output_folder,['Elastix_' num2str(idx) '_' num2str(lev) '.txt']);
                            pars = MRT_Library.ReadElastixParameters(f2r);
                            if(lev < max_depth(idx))
                                pars = MRT_Library.SetElastixParameter(pars,'InitialTransformParametersFileName',...
                                    fullfile(pwd,output_folder,['Elastix_' num2str(idx) '_' num2str(lev+1) '.txt']));
                                MRT_Library.WriteElastixParameters(f2r,pars);
                            else
                                if(idx < length(transformation_files))
                                    pars = MRT_Library.SetElastixParameter(pars,'InitialTransformParametersFileName',...
                                        fullfile(pwd,output_folder,['Elastix_' num2str(idx+1) '_' num2str('1') '.txt']));
                                    MRT_Library.WriteElastixParameters(f2r,pars);
                                end
                            end
                        end
                    end
                end
            else
                % Atomic operation
                [start_index,linked_files] = MRT_Library.ElastixParametersTree(transformation_files,start_index);
                for idx=1:length(linked_files)
                    copyfile(linked_files{idx},...
                        fullfile(output_folder,['Elastix_' num2str(series_index) '_' num2str(idx) '.txt']));
                end
            end
            if(start_index == 1 && series_index == 1)
                % At the end, make sure to propagate the right frame of
                % reference
                [~,linked_files] = MRT_Library.ElastixParametersTree(fullfile(output_folder,'Elastix_1_1.txt'),start_index);
                options2copy = {'Size','Index','Spacing','Origin','Direction'};
                first_par_file = linked_files{1};
                last_par_file = linked_files{end};
                first_par = MRT_Library.ReadElastixParameters(first_par_file);
                last_par = MRT_Library.ReadElastixParameters(last_par_file);
                for opt=1:length(options2copy)
                    val2copy = MRT_Library.GetElastixParameter(last_par,options2copy{opt});
                    first_par = MRT_Library.SetElastixParameter(first_par,options2copy{opt},val2copy);
                end
                MRT_Library.WriteElastixParameters(first_par_file,first_par);
            end
            
        end

        % helpers for Elastix transformation files
        function [start_index,linked_files] = ElastixParametersTree(transformation_files,start_index)
            lp = MRT_Library.ReadElastixParameters(transformation_files);
            linked_files = {transformation_files};
            linked_file = MRT_Library.GetElastixParameter(lp,'InitialTransformParametersFileName');
            linked_file = strrep(linked_file,'"','');
            linked_file = strrep(linked_file,' ','');
            if(strcmpi(linked_file,'NoInitialTransform') == false && exist(linked_file,'file') < 1)
                % Try copying the filepath of the previous file
                [fp,~] = fileparts(transformation_files);
                linked_file = fullfile(fp,linked_file);
                if(exist(linked_file,'file') < 1)
                    disp('Cannot locate at least one of the linked files');
                end
            end
            disp(['Recursion depth ' num2str(start_index) ' on ' linked_file])

%             linked_files(end+1) = {linked_file};
            if(strcmpi(linked_file,'NoInitialTransform') == false)
                % There is a linked parameter to copy
                [start_index,linked_file] = MRT_Library.ElastixParametersTree(linked_file,start_index+1);
                for id=1:length(linked_file)
                    linked_files(end+1) = linked_file(id);
                end
            else
                % This is the last one
            end
        end
        function WriteElastixParameters(fout,pars)
            f = fopen(fout,'wt');
            for ix=1:length(pars)
                fprintf(f,'%s%s',pars{ix},newline);
            end
            fclose(f);
        end
        function params = ReadElastixParameters(fin)
            params = {};
            f = fopen(fin,'rt');
            while(feof(f) == false)
                params(end+1) = {fgetl(f)};
            end
            fclose(f);
        end
        function [parval,ix] = GetElastixParameter(params,parameter_name)
            parval = NaN;
            for ix=1:length(params)
                if(contains(params{ix},parameter_name))
                    sp = strfind(params{ix},' ');
                    ep = strfind(params{ix},')');
                    parval = params{ix}(sp:ep-1);
                    return
                end
            end
        end
        function params = SetElastixParameter(params,parameter_name,newvalue)
            [~,ix] = MRT_Library.GetElastixParameter(params,parameter_name);
            sp = strfind(params{ix},' ');
            ep = strfind(params{ix},')');
            nv = params{ix};
            nv = [nv(1:sp-1) ' ' newvalue nv(ep:end)];
            params(ix) = {nv};
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

