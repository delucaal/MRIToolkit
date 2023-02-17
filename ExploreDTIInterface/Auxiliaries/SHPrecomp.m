%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%
%%% Alberto De Luca - alberto@isi.uu.nl $%%%
%%% Distributed under the terms of LGPLv3  %%%



classdef SHPrecomp
% Originally written from Ben Jeurissen (ben.jeurissen@uantwerpen.be)
% under the supervision of Alexander Leemans (a.leemans@umcutrecht.nl)
%
    % Perform fast operations on SH coefficients, using precomputed
    % legendre functions. Always use SHPrecomp.init(lmax) first!!
    %
    % Copyright Ben Jeurissen (ben.jeurissen@ua.ac.be)
    %
    methods (Static = true, Access = private)
        function q = precompute(lmax, num)
            el = (0:pi/(num-1):pi);
            q = zeros(num,lmax*(lmax+4)/4+1,'single');
            sign = (-1).^(0:lmax);
            for l=0:2:lmax
                i = (l*l)/4 + 1;
                q(:,i:i+l) = sign(ones(1,num),1:l+1).*legendre(l,cos(el),'sch')'.*sqrt((2*l+1)/(4*pi));
            end
            q = q';
        end
    end
    methods (Static = true, Access = public)
        function init(lmax, num)
            global SHPrecompSettings
            if ~exist('num','var'); num = 512; end
            mex_sh('init', SHPrecomp.precompute(lmax,num));
            SHPrecompSettings.lmax = lmax;
            SHPrecompSettings.num = num;
        end
        
        function [peaks, vals] = peaks(sh, init_dir)
            [peaks, vals] = mex_sh('peaks', sh, init_dir);
        end
        
        function [peaks, vals] = all_peaks(sh, min_val, max_nr, ncores, lmax, ndirs)
            if(nargin > 3)
                [peaks,vals] = SHPrecomp.all_peaks_parallel(sh, min_val, max_nr, ncores, lmax, ndirs);
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
                [peaks_i, vals_i] = mex_sh('peaks', sh_i, init_dir);
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
            pctRunOnAll(['SHPrecomp.init(' num2str(SHPrecompSettings.lmax) ', ' num2str(SHPrecompSettings.num) ')']);
                        
            parfor (i = 1:size_sh(2),ncores)
                sh_i = sh(:,ones(1,size_init_dir(2))*i);
                [peaks_i, vals_i] = mex_sh('peaks', sh_i, init_dir);
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
            [dir, val] = mex_sh('sample', sh, init_dir, max_angle, min_val, trials);
        end
        
        function val = value(sh, dir)
            val = mex_sh('eval', sh, dir);
        end
    end
end
