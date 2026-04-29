%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

classdef SH_v2
    % Alberto De Luca (a.deluca-2@umcutrecht.nl)
    % Based on SH, previously written by Ben Jeurissen (ben.jeurissen@uantwerpen.be)
    % under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
    properties (Access = public)
        lmax
        dir
        sh
        shinv
        hat
        leverage_fact
        basis % 0 (EDTI), 1 (DESCOTEAUX07/DIPY), 2 (MRTRIX)
    end

    methods (Access = public)
        function this = SH_v2(lmax, dir, basis)
            this.lmax = lmax;
            this.dir = dir;
            if(exist('basis','var') < 1)
                basis = 0;
            end
            this.basis = basis;
            this.sh = SH_v2.eval(lmax, dir, basis);
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

        % ADL
        function coef = reg_coef(this,amp,lambda)
            if ndims(amp) ~= 2; [amp, mask] = vec(amp); end;
            B = this.sh'*this.sh;
            %             coef = this.sh\amp;
            coef = inv(B+lambda*eye(size(B,1)))*(this.sh'*amp);
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
            lmax = floor(SH_v2.n2lmax(n));
            if ~iseven(lmax)
                lmax = lmax - 1;
            end
        end

        function [f, df_del, df_daz, d2f_del2, d2f_daz2, d2f_deldaz] = eval(lmax,dir,basis)
            if ~exist('basis','var'); basis = 0; end % 0=EDTI,1=Descoteaux07,2=MRtrix

            if size(dir,2) == 3
                switch basis
                    case {1,2} % Descoteaux07 / MRtrix
                        dir = dir(:,[2 1 3]);
                        if(basis == 1)
                            dir(:,1) = -dir(:,1);
                        elseif(basis == 2)
                            dir(:,3) = -dir(:,3);
                        end
                end
                dir = c2s(dir);
            end

            if size(dir,2) == 3
                dir = c2s(dir);
            end

            el = dir(:,1);   % theta
            az = dir(:,2);   % phi
            num = size(dir,1);

            f = zeros(num, SH_v2.lmax2n(lmax), 'double');

            for l = 0:2:lmax
                Plm = legendre(l, cos(el), 'sch')'; % (num x (l+1))
                norm = sqrt((2*l+1)/(4*pi));

                switch basis

                    % =========================
                    % ===== EDTI  =====
                    % =========================
                    case 0
                        loff = l*(l+1)/2 + 1;

                        for m = 0:l
                            Ylm = (-1)^m * norm * Plm(:,m+1) .* exp(1i*m*az);

                            if m == 0
                                f(:,loff) = real(Ylm);
                            else
                                f(:,loff-m) = imag(Ylm);
                                f(:,loff+m) = real(Ylm);
                            end
                        end

                    case 1
                        % =========================
                        % === DESCOTEAUX07/DIPY ===
                        % =========================
                        for m = -l:l
                            % Se Plm è 'sch', ha già dentro sqrt(2 * (l-m)!/(l+m)!)
                            % La costante di normalizzazione standard diventa semplicissima:
                            norm_const = sqrt((2*l + 1) / (4 * pi));

                            % Se m > 0, Schmidt ha un sqrt(2) extra.
                            % Per m = 0, Schmidt è uguale a Legendre standard.
                            % Quindi per m=0 non serve correggere, per m!=0 lo sqrt(2) è già lì.

                            val_Plm = Plm(:, abs(m)+1);
                            if mod(abs(m), 2) == 1
                                val_Plm = -val_Plm; % Condon-Shortley
                            end

                            index = (l^2 + l + 2)/2 + m;
                            if m == 0
                                f(:,index) = norm_const * val_Plm;
                            elseif m > 0
                                % sin(m*az) * norm * Plm_sch
                                f(:,index) = norm_const * val_Plm .* sin(m * az);
                            else % m < 0
                                % cos(m*az) * norm * Plm_sch
                                f(:,index) = norm_const * val_Plm .* cos(abs(m) * az);
                            end
                        end
                    case 2
                        % =========================
                        % === MRTRIX ===
                        % =========================
                        loff = l*(l+1)/2 + 1;
                        % Usa 'abs' per ottenere i polinomi standard, poi applichiamo noi la norma
                        P = legendre(l, cos(el));
                        for m = 0:l
                            % Fattore di normalizzazione standard (non Schmidt)
                            % Includiamo (l-m)!/(l+m)!
                            f_norm = sqrt((2*l+1)/(4*pi) * factorial(l-m)/factorial(l+m));

                            if m == 0
                                f(:, loff) = f_norm * P(1,:)';
                            else
                                % MRtrix usa sin(m*phi) per m < 0 e cos(m*phi) per m > 0
                                % Il fattore sqrt(2) serve per l'ortonormalità delle SH reali
                                common = sqrt(2) * f_norm * P(m+1,:)';
                                f(:, loff - m) = common .* sin(m * az); % m negativo
                                f(:, loff + m) = common .* cos(m * az); % m positivo
                            end
                        end
                    case 21
                        % =========================
                        % === MRTRIX ===
                        % =========================
                        loff = l*(l+1)/2 + 1;

                        for m = 0:l

                            % complex SH (NO extra (-1)^m here!)
                            Ylm = norm * Plm(:,m+1) .* exp(1i*m*az);

                            if m == 0
                                % m = 0 is purely real
                                f(:, loff) = real(Ylm);

                            else
                                % MRtrix-style real basis:
                                % symmetric combination of ±m

                                f(:, loff - m) = sqrt(2) * imag(Ylm);
                                f(:, loff + m) = sqrt(2) * real(Ylm);
                            end
                        end
                end
            end

            if nargout > 1 % Numerical derivatives instead of analytical as before
                eps = 1e-5;

                df_del = zeros(size(f));
                df_daz = zeros(size(f));

                if nargout > 3
                    d2f_del2 = zeros(size(f));
                    d2f_daz2 = zeros(size(f));
                    d2f_deldaz = zeros(size(f));
                end

                for i = 1:num
                    theta = el(i);
                    phi   = az(i);

                    % --- first derivatives ---
                    f_p = SH_v2.eval(lmax, [theta+eps phi], basis);
                    f_m = SH_v2.eval(lmax, [theta-eps phi], basis);
                    df_del(i,:) = (f_p - f_m) / (2*eps);

                    f_p = SH_v2.eval(lmax, [theta phi+eps], basis);
                    f_m = SH_v2.eval(lmax, [theta phi-eps], basis);
                    df_daz(i,:) = (f_p - f_m) / (2*eps);

                    if nargout > 3
                        % --- second derivatives ---
                        f0 = f(i,:);

                        f_p = SH_v2.eval(lmax, [theta+eps phi], basis);
                        f_m = SH_v2.eval(lmax, [theta-eps phi], basis);
                        d2f_del2(i,:) = (f_p - 2*f0 + f_m) / (eps^2);

                        f_p = SH_v2.eval(lmax, [theta phi+eps], basis);
                        f_m = SH_v2.eval(lmax, [theta phi-eps], basis);
                        d2f_daz2(i,:) = (f_p - 2*f0 + f_m) / (eps^2);

                        % mixed derivative
                        f_pp = SH_v2.eval(lmax, [theta+eps phi+eps], basis);
                        f_pm = SH_v2.eval(lmax, [theta+eps phi-eps], basis);
                        f_mp = SH_v2.eval(lmax, [theta-eps phi+eps], basis);
                        f_mm = SH_v2.eval(lmax, [theta-eps phi-eps], basis);

                        d2f_deldaz(i,:) = (f_pp - f_pm - f_mp + f_mm) / (4*eps^2);
                    end
                end
            end
        end
    end
end