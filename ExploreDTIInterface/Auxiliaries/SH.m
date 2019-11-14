%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Distributed under the terms of LGPLv3  %%%
classdef SH
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)    
    properties (Access = private)
        lmax
        dir
        sh
        shinv
        hat
        leverage_fact
    end
    
    methods (Access = public)
        function this = SH(lmax, dir)
            this.lmax = lmax;
            this.dir = dir;
            this.sh = SH.eval(lmax, dir);
            this.shinv = pinv(this.sh);
            this.hat = this.sh*this.shinv;
            if sum(real(sqrt(1-diag(this.hat))))<10^-1
                this.leverage_fact = 1;
            else
                this.leverage_fact = 1./sqrt(1-diag(this.hat));
            end
        end
        
        function coef = coef(this, amp)
            if ndims(amp) ~= 2; [amp, mask] = vec(amp); end;
            coef = this.shinv*amp;
            if exist('mask','var'); coef = unvec(coef,mask); end;
        end
        
        function amp = amp(this, coef)
            if ndims(coef) ~= 2; [coef, mask] = vec(coef); end;
            amp = this.sh*coef;
            if exist('mask','var'); amp = unvec(amp,mask); end;
        end
        
        function [amp_fit, mod_res] = fit(this, amp)
            amp_fit = this.hat*amp; % calculate SH fit
            if nargout > 1
                res = amp-amp_fit; % calculate raw residual
                mod_res = bsxfun(@times,res,this.leverage_fact); % modify residual for leverage
                mod_res = bsxfun(@minus,mod_res,mean(mod_res)); % center modified residual around mean
            end
        end
    end
    
    methods (Access = public, Static = true)
        function n = lmax2n(lmax)
            n = (lmax+1).*(lmax+2)./2;
        end
        
        function lmax = n2lmax(n)
            lmax = 2.*(floor((sqrt(1+8.*n)-3)./4));
        end
        
        function lmax = maxlmax(n)
            lmax = floor(SH.n2lmax(n));
            if ~iseven(lmax)
                lmax = lmax - 1;
            end
        end
        
        function [f, df_del, df_daz, d2f_del2, d2f_daz2, d2f_deldaz] = eval(lmax,dir)
            if size(dir,2) == 3
                dir = c2s(dir);
            end
            el = dir(:,1)';
            az = dir(:,2);
            
            num = size(dir,1);
            
            f = zeros(num,SH.lmax2n(lmax),'double');
            if nargout > 1
                df_del = zeros(size(f),'double');
                df_daz = zeros(size(f),'double');
                if nargout > 3
                    d2f_del2 = zeros(size(f),'double');
                    d2f_daz2 = zeros(size(f),'double');
                    d2f_deldaz = zeros(size(f),'double');
                end
            end
            sign = (-1).^(0:lmax);
            for l=0:2:lmax
                q = sign(ones(1,num),1:l+1);
                q = q.*legendre(l,cos(el),'sch')';
                q = q.*sqrt((2*l+1)/(4*pi));
                loff = l*(l+1)/2 + 1;
                for m=0:l
                    if m
                        cosmaz = cos(m*az); sinmaz = sin(m*az);
                        f(:,loff+m) = q(:,m+1).*cosmaz;
                        f(:,loff-m) = q(:,m+1).*sinmaz;
                        if nargout > 1
                            if m == 1
                                tmp = sqrt((l+m)*(l-m+1)*2).*q(:,m);
                            else
                                tmp = sqrt((l+m)*(l-m+1)  ).*q(:,m);
                            end
                            if m < l
                                tmp = tmp - sqrt((l-m)*(l+m+1)).*q(:,m+2);
                            end
                            tmp = -tmp/2;
                            if nargout > 3
                                tmp2 = -((l+m)*(l-m+1) + (l-m)*(l+m+1)).*q(:,m+1);
                                if (m == 1)
                                    tmp2 = tmp2 - ((l+1)*l).*q(:,2);
                                else
                                    if (m == 2)
                                        tmp2 = tmp2 + sqrt((l+m)*(l-m+1)*(l+m-1)*(l-m+2)*2).*q(:,m-1);
                                    else
                                        tmp2 = tmp2 + sqrt((l+m)*(l-m+1)*(l+m-1)*(l-m+2)).*q(:,m-1);
                                    end
                                end
                                if (l > m+1)
                                    tmp2 = tmp2 + sqrt((l-m)*(l+m+1)*(l-m-1)*(l+m+2)).*q(:,m+3);
                                end
                                tmp2 = tmp2/4.0;
                                
                                d2f_del2(:,loff+m) = tmp2.*cosmaz;
                                d2f_del2(:,loff-m) = tmp2.*sinmaz;
                                d2f_daz2(:,loff+m) = -q(:,m+1).*cosmaz.*m.*m;
                                d2f_daz2(:,loff-m) = -q(:,m+1).*sinmaz.*m.*m;
                                d2f_deldaz(:,loff+m) = -tmp.*sinmaz.*m;
                                d2f_deldaz(:,loff-m) =  tmp.*cosmaz.*m;
                            end
                            df_del(:,loff+m) = tmp.*cosmaz;
                            df_del(:,loff-m) = tmp.*sinmaz;
                            
                            df_daz(:,loff+m) = -q(:,m+1).*sinmaz.*m;
                            df_daz(:,loff-m) =  q(:,m+1).*cosmaz.*m;
                        end
                    else
                        f(:,loff) = q(:,1);
                        if l > 1
                            if nargout > 1
                                df_del(:,loff) = q(:,2).*sqrt(l*(l+1)/2);
                                if nargout > 3
                                    d2f_del2(:,loff) = (sqrt(l*(l+1)*(l-1)*(l+2)/2) * q(:,3) - l*(l+1) * q(:,1))/2;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end