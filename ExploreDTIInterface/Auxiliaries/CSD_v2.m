%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

classdef CSD_v2
% Alberto De Luca a.deluca-2@umcutrecht.nl
% Based on code previosuly written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
%
    % An updated version of the Constrained Spherical Deconvolution (CSD) class
    % This version supports multiple SH definitions
    % Usage:
    % r_sh = CSD_v2.response(dwi, grad, bval, minFA, lmax);
    % csd = CSD_v2(r_sh, target_lmax, basis); % 0 = EDTI, 1 = DIPY, 2 =
    % MRTRIX
    % fod_sh = csd.deconv(dwi_sh);
    %
    % for super-resolution: make sure target_lmax is higher than SH order of r_sh
    %
    % default parameters can be changed after construction by:
    % csd.setLambda(lambda);
    % csd.setTau(tau)
    % csd.setIter(iter)
    % csd.setInitLmax(init_lmax)
    % csd.setNegativeDirs(dirs)
    %
    %
    properties (GetAccess = public, SetAccess = private)
        lmax;
        basis; %0 EDTI, 1 DIPY, 2 MRTRIX
    end
    
    properties (Access = private)
        n;
        convmat;
        hr_sh;
        lambda_;
        r_rh;
        my_fod;
        
        lambda = 1;
        tau = 0.1;
        dirs = textread('dir300.txt');
        niter = 50;
        init_lmax = 4;
    end
    
    methods (Access = public)
        function this = CSD_v2(r_sh, target_lmax, basis)
            if ~exist('target_lmax','var') || isempty(target_lmax)
                target_lmax = SH_v2.n2lmax(size(r_sh,1));
            end
            if(exist('basis','var') < 1)
                basis = 0; % Default to EDTI
            end
            this.basis = basis;
            this.lmax = target_lmax;
            r_sh_n = size(r_sh,1);
            r_sh_lmax = SH_v2.n2lmax(r_sh_n);
            target_n = SH_v2.lmax2n(target_lmax);
            
            lmax_ = min(r_sh_lmax,target_lmax);
            this.n = min(r_sh_n,target_n);
            
            d_sh = SH_v2.eval(lmax_,[0 0],this.basis)';
            k = find(d_sh);
            this.r_rh = r_sh(k)./d_sh(k);
            m = [];
            for l = 0:2:lmax_
                m = [m; ones(2*l+1,1)*this.r_rh(l/2+1)];
            end
            this.convmat = diag(m);
            
            this.hr_sh = SH_v2.eval(this.lmax,c2s(this.dirs),this.basis);
            this.convmat = [this.convmat zeros(size(this.convmat,1),size(this.hr_sh,2)-size(this.convmat,2)) ];
            this.lambda_ = this.lambda*size(this.convmat,1)*this.r_rh(1)/size(this.hr_sh,1);
        end
        
        function this = setLambda(this, lambda)
            this.lambda = lambda;
            this.lambda_ = this.lambda*size(this.convmat,1)*this.r_rh(1)/size(this.hr_sh,1);
        end
        
        function this = setTau(this, tau)
            this.tau = tau;
        end
        
        function this = setIter(this, iter)
            this.niter = iter;
        end
        
        function this = setInitLmax(this, init_lmax)
            this.init_lmax = init_lmax;
        end
        
        function this = setNegativeDirs(this, dirs)
            this.dirs = dirs;
            this.hr_sh = SH_v2.eval(this.lmax,c2s(this.dirs),this.basis);
            this.lambda_ = this.lambda*size(this.convmat,1)*this.r_rh(1)/size(this.hr_sh,1);
        end
        
        function fod_sh = process(this, dwi_sh)
            fod_sh = this.deconv(dwi_sh);
        end
        
        function dwi_sh = conv(this, fod_sh)
            if ndims(fod_sh) ~= 2; [fod_sh, mask] = vec(fod_sh); end;
            dwi_sh = this.convmat*fod_sh;
            if exist('mask','var'); dwi_sh = unvec(dwi_sh,mask); end;
        end
        
        function fod_sh = deconv(this, dwi_sh)
            if ndims(dwi_sh) ~= 2; [dwi_sh, mask] = vec(dwi_sh); end;
            fod_sh = this.convmat\dwi_sh(1:this.n,:); 
            fod_sh(SH_v2.lmax2n(this.init_lmax)+1:end,:) = 0;
            threshold = this.tau*mean(this.hr_sh*fod_sh,1);
            fod_sh = CSD_v2.actual_deconv(fod_sh, dwi_sh(1:this.n,:), this.convmat, this.hr_sh, threshold, this.lambda_, this.niter);
            if exist('mask','var'); fod_sh = unvec(fod_sh,mask); end;
        end
    end
    
    methods (Access = private, Static = true)
        function f_sh = actual_deconv(f_sh,dwi_sh,fconv,hr_sh,threshold,lambda,niter)
            parfor i = 1:size(f_sh,2) % for each voxel
                f_sh_i = f_sh(:,i); 
                dwi_sh_i = dwi_sh(:,i); threshold_i = threshold(i);
                for it = 1:niter % maximum of 50 CSD iterations
                    f_hr = hr_sh*f_sh_i; % calculate high-res FOD
                    neg = find (f_hr < threshold_i); % find negative (and small) directions in high-res FOD
                    if size(neg,1) + size(fconv,1) < size(hr_sh,2) % if there aren't enough negative directions to allow super resolution
                        disp('Warning: too few negative directions identified - failed to converge'); % stop
                        break;
                    end
                    m = [fconv; lambda*hr_sh(neg,:)]; % assume the negative directions
                    s = [dwi_sh_i; zeros(size(neg,1),1)]; % are zero
                    f_sh_i = m\s; % re-estimate f_sh as if more directions are available
                    ssd = sum((f_sh_i-f_sh(:,i)).^2);
                    f_sh(:,i) = f_sh_i;
                    if (ssd == 0.0) % if there is no improvement of f_sh over the previous iteration
                        break; % we have converged
                    end
                end
            end
        end
    end
    
    methods (Access = public, Static = true)
        function r_sh = response(dwi, grad, bval, minFA, lmax, basis)
            if ndims(dwi) ~= 2; dwi = vec(dwi); end;
            if(exist('basis','var') < 1)
                basis = 0;
            end

            % calculate response function
            shell = abs(grad(:,4)-bval)<1;
            shell_grad = grad(shell,1:3);
            
            %% check grad
            if size(grad,2) ~= 4
                error('grad must contain 4 columns');
            end
%             grad = single(grad);
            
            %% check minFA
            if ~exist('minFA','var') || isempty(minFA)
                minFA = 0.7;
            else
                if (minFA > 1 || minFA < 0)
                    error('minFA outside [0 1]');
                end
            end
%             minFA = single(minFA);
            %% check lmax
            lmax_grad = SH_v2.maxlmax(size(shell_grad, 1));
            %disp(['highest possible lmax is ' num2str(lmax_grad)]);
            if ~exist('lmax','var') || isempty(lmax)
                lmax = lmax_grad;
            else
                %disp(['requested lmax is ' num2str(lmax)]);
                if (lmax < 0 || lmax > lmax_grad)
                    error(['lmax should be in [0 ' num2str(lmax_grad) '] for this dataset']);
                end
                if (~iseven(lmax))
                    error('lmax should be even');
                end
            end
            %disp(['effective lmax is ' num2str(lmax)]);
            
            %% dti
            dti = DTI(grad);
            dt = dti.dt(dwi,true);
            [eigval, fe, se, te] = DTI.eig(dt);
            fa = DTI.fa(eigval);
            sf_mask = fa > minFA;
            dwi = dwi(shell,sf_mask);
            fe = fe(:,sf_mask); se = se(:,sf_mask); te = te(:,sf_mask);
            %% estimate sh coefficients of response function
            r_sh = zeros(SH_v2.lmax2n(lmax),1);
            for i = 1:size(dwi,2)
                rot_grad = zeros(size(shell_grad,1),3);
                rot = [te(:,i) se(:,i) fe(:,i)];
                for j = 1:size(shell_grad,1)
                    rot_grad(j,:) = shell_grad(j,:)*rot;
                end
                r_sh = r_sh + (SH_v2.eval(lmax,c2s(rot_grad),basis)\dwi(:,i));
            end
            r_sh = r_sh ./ size(dwi,2);
            i = 0;
            for l_ = 0:2:lmax
                for m = -l_:l_
                    i = i + 1;
                    if m ~= 0
                        r_sh(i) = 0;
                    end
                end
            end
            %disp(r_sh(r_sh~=0)')
        end        
    end
end
