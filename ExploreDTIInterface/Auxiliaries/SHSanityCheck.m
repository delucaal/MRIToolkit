%% Sanity check: SH vs SH_v2 (EDTI basis only)

clear; clc;

lmax = 8;

% random directions in Cartesian
N = 500;
dir = randn(N,3);
dir = dir ./ vecnorm(dir,2,2);

% convert to single (important: both pipelines use single in MEX)
dir_single = single(dir');

% fake SH coefficients (consistent size)
ncoef = SH.lmax2n(lmax);
sh = single(randn(ncoef,N));

%% --- Evaluate old MATLAB SH ---
basis = 0; % EDTI

tic
[f_old, df_del_old, df_daz_old] = SH.eval(lmax, dir);
toc

%% --- Evaluate new MATLAB SH ---
basis = 0; % EDTI

tic
[f_new, df_del_new, df_daz_new] = SH_v2.eval(lmax, dir, basis);
toc

%% Compare

figure
subplot(311)
plot(f_old,'*b')
hold on
plot(f_new,'or')
subplot(312)
plot(df_del_old,'*b')
hold on
plot(df_del_new,'or')
subplot(313)
plot(df_daz_old,'*b')
hold on
plot(df_daz_new,'or')

disp(['Maximum error on SH: ' num2str(max(abs(f_new(:) - f_old(:))))])
disp(['Maximum error on dSH/del: ' num2str(max(abs(df_del_new(:) - df_del_old(:))))])
disp(['Maximum error on dSH/daz: ' num2str(max(abs(df_daz_new(:) - df_daz_old(:))))])

%% --- Evaluate new MEX SH_v2 ---
lmax = 8;
num = 512;
for basis = 0:2

    SHPrecomp_v2.init(lmax,num,basis);
    
    sh = SH_v2(lmax,dir,basis);
    
    c = single(randn(45,1));               % coefficients
    
    f1 = sh.amp(c);                   % MATLAB reference
    c_mat = repmat(c, 1, size(dir,1));   % 45 × 512
    f2 = mex_sh_v2('eval', c_mat, single(dir')); % MEX
    
    max(abs(f1 - f2'))
    
    figure
    plot(f1,'*b')
    hold on
    plot(f2,'or')
end

