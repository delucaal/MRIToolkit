%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



classdef SD
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
%
    % Spherical Deconvolution (SD) class
    %
    % Usage:
    % r_sh = SD.response(dwi, grad, bval, minFA, lmax);
    % sd = SD(r_sh, target_lmax);
    % fod_sh = sd.deconv(dwi_sh);
    %
    % see also CSD
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    properties (GetAccess = public, SetAccess = private)
        lmax;
    end
    
    properties (Access = private)
        convmat;
        n;
    end
    
    methods (Access = public)
        function this = SD(r_sh, target_lmax)
            r_sh_n = size(r_sh,1);
            this.lmax = target_lmax;
            this.n = SH.lmax2n(this.lmax);
            if (this.n > r_sh_n); disp('Error: SH order of response function not high enough for target lmax'); error(''); end;
            d_sh = SH.eval(this.lmax,[0 0])';
            k = find(d_sh);
            r_rh = r_sh(k)./d_sh(k);
            m = [];
            for l = 0:2:this.lmax
                m = [m; ones(2*l+1,1)*r_rh(l/2+1)];
            end
            this.convmat = diag(m);
        end
        
        function fod_sh = process(this, dwi_sh)
            fod_sh = this.deconv(dwi_sh);
        end
        
        function fod_sh = deconv(this, dwi_sh)
            if ndims(dwi_sh) ~= 2; [dwi_sh, mask] = vec(dwi_sh); end;
            if size(dwi_sh,1) < this.n; error('SH order of dwi not high enough for target lmax'); end;
            fod_sh = this.convmat\dwi_sh(1:this.n,:);
            if exist('mask','var'); fod_sh = unvec(fod_sh,mask); end;
        end
    end
        
    methods (Access = public, Static = true)
        function r_sh = response(dwi, grad, bval, minFA, lmax)
            if ndims(dwi) ~= 2; dwi = vec(dwi); end;
            
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
            lmax_grad = SH.maxlmax(size(shell_grad, 1));
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
            r_sh = zeros(SH.lmax2n(lmax),1);
            for i = 1:size(dwi,2)
                rot_grad = zeros(size(shell_grad,1),3);
                rot = [te(:,i) se(:,i) fe(:,i)];
                for j = 1:size(shell_grad,1)
                    rot_grad(j,:) = shell_grad(j,:)*rot;
                end
                r_sh = r_sh + (SH.eval(lmax,c2s(rot_grad))\dwi(:,i));
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