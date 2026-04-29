%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%

classdef SHPrecomp_v2
    % Alberto De Luca (a.deluca-2@umcutrecht.nl)
    % Based on SHPrecomp_v2 from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
    % under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
    %
    % Perform fast operations on SH coefficients, using precomputed
    % legendre functions. Always use SHPrecomp_v2.init(lmax) first!!
    % This version supports multiple SH definitions (EDTI, MRTRIX, DIPY)
    %
    methods (Static = true, Access = private)
        % function q = precompute(lmax, num, basis)
        %     if ~exist('basis','var'); basis = 0; end
        %     el = (0:pi/(num-1):pi);
        %     q = zeros(num,lmax*(lmax+4)/4+1,'single');
        %     if(basis == 0 || basis == 1)
        %         sign = (-1).^(0:lmax);
        %     else
        %         % MRTRIX
        %         sign = 1.^(0:lmax);
        %     end
        %     for l=0:2:lmax
        %         i = (l*l)/4 + 1;
        %         q(:,i:i+l) = sign(ones(1,num),1:l+1).*legendre(l,cos(el),'sch')'.*sqrt((2*l+1)/(4*pi));
        %     end
        %     q = q';
        % end
        function q = precompute(lmax, num, basis)
            if ~exist('basis','var'); basis = 0; end

            % Generiamo i punti di elevazione (theta)
            el = (0:pi/(num-1):pi);

            % Numero di coefficienti nella tabella precomputata (solo m >= 0)
            num_coeffs = (lmax/2 + 1)^2; % Equivalente a lmax*(lmax+4)/4 + 1
            q = zeros(num, num_coeffs, 'single');

            for l = 0:2:lmax
                % Indice di inizio per il grado l nella tabella m-positiva
                idx_start = (l*l)/4 + 1;

                % Calcoliamo i polinomi di Legendre associati
                % Usiamo 'sch' per semplicità, ma con correzioni per MRtrix

                switch basis
                    case {0, 1} % EDTI o Descoteaux/DIPY
                        Plm = legendre(l, cos(el), 'sch')';
        
                        % Normalizzazione standard comune
                        norm_const = sqrt((2*l+1)/(4*pi));
                        % Applichiamo la fase di Condon-Shortley (-1)^m
                        m_indices = 0:l;
                        s = (-1).^m_indices;
                        q(:, idx_start:idx_start+l) = Plm .* (s .* norm_const);

                    case 2 % MRTRIX
                        % MRtrix NON vuole la fase (-1)^m
                        % Ma vuole sqrt(2) per tutti i termini m > 0
                        P = legendre(l, cos(el));
                        for m = 0:l
                            f_norm = sqrt(((2*l+1)/(4*pi)) * (factorial(l-m)/factorial(l+m)));
                            val = P(m+1, :) * f_norm;
                            if m > 0
                                val = val * sqrt(2);
                            end
                            q(:, idx_start + m) = single(val');
                        end
                end
            end

            % Trasponiamo per avere [Coefficienti x Direzioni] come atteso dal MEX
            q = q';
        end

    end
    methods (Static = true, Access = public)
        function init(lmax, num, basis)
            global SHPrecompSettings
            if ~exist('num','var'); num = 512; end
            if ~exist('basis','var'); basis = 0; end
            mex_sh_v2('init', SHPrecomp_v2.precompute(lmax,num,basis),basis);
            SHPrecompSettings.lmax = lmax;
            SHPrecompSettings.num = num;
            SHPrecompSettings.basis = basis;
        end

        function [peaks, vals] = peaks(sh, init_dir)
            [peaks, vals] = mex_sh_v2('peaks', sh, init_dir);
        end

        function [peaks, vals] = all_peaks(sh, min_val, max_nr, ncores, lmax, ndirs)
            if(nargin > 3)
                [peaks,vals] = SHPrecomp_v2.all_peaks_parallel(sh, min_val, max_nr, ncores, lmax, ndirs);
                return
            end

            if ndims(sh) ~= 2; [sh, mask] = vec(sh); end;
            %             init_dir = single(textread('jones30.txt')');
            init_dir = single(textread('dir060.txt')');
            sh = single(sh);
            size_sh = size(sh);
            size_init_dir = size(init_dir);

            peaks = cell(1,size_sh(2));
            vals = cell(1,size_sh(2));

            for i=1:size_sh(2)
                sh_i = sh(:,ones(1,size_init_dir(2))*i);
                [peaks_i, vals_i] = mex_sh_v2('peaks', sh_i, init_dir);
                mask_t = vals_i > min_val; peaks_i = peaks_i(:,mask_t); vals_i = vals_i(:,mask_t);

                if size(peaks_i,2) > 0
                    peaks_i_ok = NaN(3,max_nr,'single');
                    vals_i_ok = NaN(1,max_nr,'single');
                    [vals_i,index] = sort(vals_i, 'descend');
                    peaks_i = peaks_i(:,index);

                    pcount = 1;
                    peaks_i_ok(:,pcount) = peaks_i(:,1);
                    vals_i_ok(:,pcount) = vals_i(:,1);

                    for j = 2:size(peaks_i,2)
                        if pcount == max_nr
                            break;
                        end
                        ok = true;
                        for k = 1:pcount
                            dt = sum(peaks_i_ok(1:3,k).*peaks_i(:,j),1);
                            angle = (180/pi)*real(acos(abs(dt)));
                            if angle < 10
                                ok = false;
                                break;
                            end
                        end
                        if ok
                            pcount = pcount + 1;
                            peaks_i_ok(:,pcount) = peaks_i(:,j);
                            vals_i_ok(:,pcount) = vals_i(:,j);
                        end
                    end
                    peaks{i} = peaks_i_ok(:,1:pcount);
                    vals{i} = vals_i_ok(:,1:pcount);
                end
            end
            if exist('mask','var');
                peaks_ = cell(size(mask)); peaks_(mask) = peaks; peaks = peaks_;
                vals_ = cell(size(mask)); vals_(mask) = vals; vals = vals_;
            end;
        end

        function [peaks, vals] = all_peaks_parallel(sh, min_val, max_nr, ncores, lmax, ndirs)
            global SHPrecompSettings

            if ndims(sh) ~= 2; [sh, mask] = vec(sh); end;
            %             init_dir = single(textread('jones30.txt')');
            init_dir = single(textread('dir060.txt')');
            sh = single(sh);
            size_sh = size(sh);
            size_init_dir = size(init_dir);

            peaks = cell(1,size_sh(2));
            vals = cell(1,size_sh(2));

            if(exist('lmax','var') < 1)
                SHPrecompSettings.lmax = SH.n2lmax(size_sh(1));
            else
                SHPrecompSettings.lmax = lmax;
            end
            SHPrecompSettings.num = ndirs;
            if(isempty(gcp('nocreate')))
                parpool(ncores)
            end
            pctRunOnAll(['SHPrecomp_v2.init(' num2str(SHPrecompSettings.lmax) ', ' num2str(SHPrecompSettings.num) ',' num2str(SHPrecompSettings.basis) ');']);

            parfor (i = 1:size_sh(2),ncores)
                sh_i = sh(:,ones(1,size_init_dir(2))*i);
                [peaks_i, vals_i] = mex_sh_v2('peaks', sh_i, init_dir);
                mask_t = vals_i > min_val; peaks_i = peaks_i(:,mask_t); vals_i = vals_i(:,mask_t);

                if size(peaks_i,2) > 0
                    peaks_i_ok = NaN(3,max_nr,'single');
                    vals_i_ok = NaN(1,max_nr,'single');
                    [vals_i,index] = sort(vals_i, 'descend');
                    peaks_i = peaks_i(:,index);

                    pcount = 1;
                    peaks_i_ok(:,pcount) = peaks_i(:,1);
                    vals_i_ok(:,pcount) = vals_i(:,1);

                    for j = 2:size(peaks_i,2)
                        if pcount == max_nr
                            break;
                        end
                        ok = true;
                        for k = 1:pcount
                            dt = sum(peaks_i_ok(1:3,k).*peaks_i(:,j),1);
                            angle = (180/pi)*real(acos(abs(dt)));
                            if angle < 10
                                ok = false;
                                break;
                            end
                        end
                        if ok
                            pcount = pcount + 1;
                            peaks_i_ok(:,pcount) = peaks_i(:,j);
                            vals_i_ok(:,pcount) = vals_i(:,j);
                        end
                    end
                    peaks{i} = peaks_i_ok(:,1:pcount);
                    vals{i} = vals_i_ok(:,1:pcount);
                end
            end
            if exist('mask','var');
                peaks_ = cell(size(mask)); peaks_(mask) = peaks; peaks = peaks_;
                vals_ = cell(size(mask)); vals_(mask) = vals; vals = vals_;
            end;
        end

        function [dir, val] = sample(sh, init_dir, max_angle, min_val, trials)
            [dir, val] = mex_sh_v2('sample', sh, init_dir, max_angle, min_val, trials);
        end

        function val = value(sh, dir)
            val = mex_sh_v2('eval', sh, dir);
        end
    end
end
