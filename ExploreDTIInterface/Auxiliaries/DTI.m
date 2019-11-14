%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
classdef DTI
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
%
	% Diffusion Tensor Imaging (DTI) class
    %
    % Usage:
    % dti = DTI(grad);
    % [dt, b0] = dti.dt(dwi);
    % [eigval, fe, se, te] = DTI.eig(dt);
    % md = DTI.md(eigval);
    % fa = DTI.fa(eigval);
    % [cl, cp, cs] = DTI.westin(eigval)
    % fefa = DTI.fefa(fe, fa);
    % adc = DTI.adc(dt, dir)
    % dwi = dti.dwi(dt, b0)
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    properties (Constant = true)
        ind = [ 1     1
            1     2
            1     3
            2     2
            2     3
            3     3];
        cnt = [ 1     2     2     1     2     1 ];
    end
    
    properties (GetAccess = public, SetAccess = private)
        b
        grad
    end
    
    methods (Access = public)
        function this = DTI(grad)
            this.grad = single(grad);
            this.b = [ ones(size(this.grad,1),1) -(this.grad(:,ones(1,6)*4).*this.grad(:,DTI.ind(1:6,1)).*this.grad(:,DTI.ind(1:6,2)))*diag(DTI.cnt)];
        end
        
        function [dt, b0] = dt(this, dwi, wls)
            dwi = single(dwi);
            if ndims(dwi) ~= 2; [dwi, mask] = vec(dwi); end;
            if nargin < 3
                wls = true;
            end
            dwi(dwi<1)=1;
            logdwi = log(dwi);
            dt = this.b\logdwi;
            if wls
                dt = mex_wls(this.b,logdwi,exp(this.b*dt));
            end
            b0 = exp(min(15,dt(1,:)));
            dt = dt(2:end,:);
            if exist('mask','var'); dt = unvec(dt,mask); b0 = unvec(b0,mask); end;
        end
        
        function [log_dwi_fit, mod_res, w] = fit(this, dwi)
            dwi = single(dwi);
            if ndims(dwi) ~= 2; [dwi, mask] = vec(dwi); end;
            dwi(dwi<1)=1;
            log_dwi = log(dwi);
            [dt, w, h] = mex_dti_wls_fit(log_dwi,this.b);
            log_dwi_fit = this.b*dt;
            mod_res = w.*(log_dwi-log_dwi_fit)./sqrt(1-h);
            mod_res = bsxfun(@minus,mod_res,mean(mod_res));
            if exist('mask','var'); log_dwi_fit = unvec(log_dwi_fit,mask); mod_res = unvec(mod_res,mask); w = unvec(w,mask); end;
        end

        function dwi = dwi(this, dt, b0)
            if ndims(dt) ~= 2; [dt, mask] = vec(dt); b0 = vec(b0); end;
            dwi = exp(this.b*[log(b0); dt]);
            if exist('mask','var'); dwi = unvec(dwi,mask); end;
        end
        
        function dwi = sim_dwi(this, fa, diffusivity, el, az, b0)
            dwi = this.dwi(DTI.sim_dt(fa, diffusivity, el, az), repmat(b0, [1 size(el,1)]));
        end
    end
    
    methods (Access = public, Static = true)
        function [eigval, fe, se, te] = eig(dt)
            if ndims(dt) ~= 2; [dt, mask] = vec(dt); end;
            dt = dt(1:6,:);
            [eigval, eigvect] = mex_dti_eig(dt);
            te = eigvect(1:3,:); se = eigvect(4:6,:); fe = eigvect(7:9,:);
            eigval = abs(eigval([3 2 1],:));
            if exist('mask','var'); eigval = unvec(eigval, mask); fe = unvec(fe, mask); se = unvec(se, mask); te = unvec(te, mask); end;
        end
        
        function fa = fa(eigval)
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            l1 = eigval(1,:); l2 = eigval(2,:); l3 = eigval(3,:);
            fa = sqrt(1/2).*sqrt((l1-l2).^2+(l2-l3).^2+(l3-l1).^2)./sqrt(l1.^2+l2.^2+l3.^2);
            if exist('mask','var'); fa = unvec(fa, mask); end;
        end
        
        function md = md(eigval)
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            md = sum(eigval,1)./3;
            if exist('mask','var'); md = unvec(md, mask); end;
        end
        
        function fefa = fefa(fe, fa)
            if ndims(fe) ~= 2; [fe, mask] = vec(fe); fa = vec(fa); end;
            fefa = abs(fe([2 1 3],:)).*fa([1 1 1],:);
            if exist('mask','var'); fefa = unvec(fefa, mask); end;
        end
        
        function adc = adc(dt, dir)
            if ndims(dt) ~= 2; [dt, mask] = vec(dt); end;
            adc = (dir(:,DTI.ind(1:6,1)).*dir(:,DTI.ind(1:6,2))) * diag(DTI.cnt) * dt;
            if exist('mask','var'); adc = unvec(adc, mask); end;
        end
        
        function [cl, cp, cs] = westin(eigval)
            if ndims(eigval) ~= 2; [eigval, mask] = vec(eigval); end;
            l1 = eigval(1,:); l2 = eigval(2,:); l3 = eigval(3,:);
            cl = (l1-l2)./l1; cp = (l2-l3)./l1; cs = l3./l1;
            if exist('mask','var'); cl = unvec(cl, mask); cp = unvec(cp, mask); cs = unvec(cs, mask); end;
        end
        
%         function odf = odf(dt, dir, k)
%             if ndims(dt) ~= 2; [dt, mask] = vec(dt); end;
%             dt = dt(1:6,:);
%             points = equatorPoints(dir,k);
%             odf = zeros(size(dir,1),size(dt,2),class(dt));
%             for i = 1:size(dir,1)
%                 p_i = squeeze(points(i,:,:));
%                 adc = DTI.adc(dt,p_i);
%                 odf(i,:) = sum(3./adc,1);
%             end
%             if exist('mask','var'); odf = unvec(odf, mask); end;
%         end
        
        function dt = sim_dt(fa, diffusivity, el, az)
            a = fa/sqrt(3-2*fa^2);
            dt_ = diffusivity*[ 1-a 0 0; 0 1-a 0; 0 0 1+2*a ];

            s = size(el);
            sinaz = sin(az);
            cosaz = cos(az);
            sinel = sin(el);
            cosel = cos(el);
            zeross = zeros(s);
            oness = ones(s);
            
            R_az = reshape([cosaz -sinaz zeross; sinaz cosaz zeross; zeross zeross oness],[s(1) 3 3]);
            R_el = reshape([cosel zeross sinel; zeross oness zeross; -sinel zeross cosel],[s(1) 3 3]);
            
            dt = zeros(6, s(1));
            for i = 1:s(1)
                R_az_i = squeeze(R_az(i,:,:));
                R_el_i = squeeze(R_el(i,:,:));
                dt_i = R_az_i*R_el_i*dt_*R_el_i'*R_az_i';
                dt(:,i) = dt_i([1 2 3 5 6 9]);
            end
        end
        
        function dt = sim_dt_eigvals(eigvals, el, az)
            dt_ = diag(eigvals);

            s = size(el);
            sinaz = sin(az);
            cosaz = cos(az);
            sinel = sin(el);
            cosel = cos(el);
            zeross = zeros(s);
            oness = ones(s);
            
            R_az = reshape([cosaz -sinaz zeross; sinaz cosaz zeross; zeross zeross oness],[s(1) 3 3]);
            R_el = reshape([cosel zeross sinel; zeross oness zeross; -sinel zeross cosel],[s(1) 3 3]);
            
            dt = zeros(6, s(1));
            for i = 1:s(1)
                R_az_i = squeeze(R_az(i,:,:));
                R_el_i = squeeze(R_el(i,:,:));
                dt_i = R_az_i*R_el_i*dt_*R_el_i'*R_az_i';
                dt(:,i) = dt_i([1 2 3 5 6 9]);
            end
        end
        
    end
end