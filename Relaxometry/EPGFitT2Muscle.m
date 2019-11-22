%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%



% A. De Luca - UMC Utrecht - alberto@isi.uu.nl
% This function performs an EPG fit of multi-echo CPMG data and is optimized
% for skeletal muscle use. The EPG fit is performed using an exhaustive grid-search method.
% The signal is represented with a bi-exponential model, assuming it to
% originate from a mixture of water and fat. The fat muscle is first fit in
% an heuristic mask of the subcutaneous fat and fixed for the voxel-wise
% fit.
% Input:
% T2MultiEchoNifti: The file path to a 4D Nifti file containing the CPMG
% series, e.g. 'T2MultiEcho.nii'
% TE: The vector of the acquired echo-times, e.g. [0.05 0.07 0.1]... (in
% seconds)
% This code uses an EPG implementation from Kelvin Layton (klayton@unimelb.edu.au)
% See below for credits
function epg_fit = EPGFitT2Muscle(T2MultiEchoNifti, TE)

disp(['Working on ' T2MultiEchoNifti]);
T2Series = T2MultiEchoNifti;

try
    T2SeriesData = load_untouch_nii(T2Series);
catch
    error('Error loading the T2 NIFTI');
end

% Discard spurious saved points
GoodPoints = TE~=0;
TE = TE(GoodPoints);

% As column vector
if(size(TE,1) < size(TE,2))
    TE = TE';
end

flip_angle_grid = 160:5:200; % The effective flip-angle (making it large results in unstable solutions)
nEchoes = length(TE);
T1w = 1.400; % s (The hypothesized T1 of water)
T1f = 0.365; % s (The hypothesized T1 of fat)
R1w = 1./T1w;
R1f = 1./T1f;

% 1st step solution grid (for subcutaneous fat)
T2w_grid = (10:5:60)*1e-3;
T2f_grid = (110:10:500)*1e-3;

% 2nd step solution grid (voxel-wise, for water)
T2w_finer_grid = (10:0.1:60)*1e-3; % The final solution of T2 water is determined in this grid

tau = double(TE(2)-TE(1));%%10e-3; % deltaTE

the_water_dictionary = zeros(nEchoes,length(T2w_grid),length(flip_angle_grid));
the_finer_water_dictionary = zeros(nEchoes,length(T2w_finer_grid),length(flip_angle_grid));
the_fat_dictionary = zeros(nEchoes,length(T2f_grid),length(flip_angle_grid));

% BUILD DICTIONARIES

for flip_i=1:length(flip_angle_grid)
    flip_angle = flip_angle_grid(flip_i)*pi/180;
    for t2w_i=1:length(T2w_grid)
        t2w = T2w_grid(t2w_i);
        sw = epg(nEchoes,tau,R1w,1/t2w,flip_angle);
        the_water_dictionary(:,t2w_i,flip_i) = sw;
    end
    for t2f_i=1:length(T2f_grid)
        t2f = T2f_grid(t2f_i);
        sf = epg(nEchoes,tau,R1f,1/t2f,flip_angle);
        the_fat_dictionary(:,t2f_i,flip_i) = sf;
    end
end

for flip_i=1:length(flip_angle_grid)
    flip_angle = flip_angle_grid(flip_i)*pi/180;
    for t2w_i=1:length(T2w_finer_grid)
        t2w = T2w_finer_grid(t2w_i);
        sw = epg(nEchoes,tau,R1w,1/t2w,flip_angle);
        the_finer_water_dictionary(:,t2w_i,flip_i) = sw;
    end
end

% Determine voxels without spurious points
Points2Fit=find(sum(T2SeriesData.img==0,4) < 2);
numOfVoxels = length(Points2Fit);

current_split_start = 1;

fw_map = zeros(size(Points2Fit));
ff_map = zeros(size(Points2Fit));
flip_angle_map = zeros(size(Points2Fit));
f_map = zeros(size(Points2Fit));
T2w_map = zeros(size(Points2Fit));
T2f_map = zeros(size(Points2Fit));
COST = zeros(size(Points2Fit));

% Calibrate fat in SC - this is done on the average SC signal
SCFatMask = GetSubcutaneousFat(T2SeriesData,length(TE)); % Compute an heuristic mask

distances = nan(length(T2w_grid),length(T2f_grid),length(flip_angle_grid));
fractions = nan(length(T2w_grid),length(T2f_grid),length(flip_angle_grid));
[sx,sy,sz,st] = size(T2SeriesData.img);
t2temp = reshape(T2SeriesData.img,sx*sy*sz,st);
AvgFatSignal = double(mean(t2temp(SCFatMask(:) > 0,GoodPoints))');
clear t2temp
signal = AvgFatSignal;

for flip_i=1:length(flip_angle_grid)
    for t2w_i=1:length(T2w_grid)
        water_signal = the_water_dictionary(:,t2w_i,flip_i);
        for t2f_i=1:length(T2f_grid)
            fat_signal = the_fat_dictionary(:,t2f_i,flip_i);
            X = [water_signal fat_signal];
            fs = X\signal;
            idx = sum(fs<0);
            distances(t2w_i,t2f_i,flip_i) = sum((signal-X*fs).^2)+idx*1e6; % Penalize negative fractions
            fractions(t2w_i,t2f_i,flip_i) = fs(1)/sum(fs);
        end
    end
end

[COST,IX] = min(distances(:));
[xi,yi,zi] = ind2sub(size(distances),IX);
FinalFatSignal = X*fs/sum(fs);

sc_fat.T2w = T2w_grid(xi);
sc_fat.T21 = T2f_grid(yi);
sc_fat.f = fs(2)/sum(fs);
sc_fat.avg_fa = flip_angle_grid(zi);
sc_fat.signal = FinalFatSignal;
disp(['SC fat T2: ' num2str(sc_fat.T21)]);
% End calibrate fat

[sx,sy,sz,st] = size(T2SeriesData.img);
T2SeriesData.img = double(T2SeriesData.img/max(T2SeriesData.img(:)));
T2SeriesData.img = reshape(T2SeriesData.img,sx*sy*sz,st);

% Voxel-wise fit
tic
parfor parind=current_split_start:numOfVoxels
    lin_ind = Points2Fit(parind);
    signal = (T2SeriesData.img(lin_ind,GoodPoints))';
    distances = nan(length(T2w_grid),length(flip_angle_grid));
    fractions = nan(length(T2w_grid),length(flip_angle_grid));
    
    if(mod(parind,10000) == 0)
        disp('Alive');
    end
    
    % Determine T2water and effective flip angle
    for flip_i=1:length(flip_angle_grid)
        for t2w_i=1:length(T2w_grid)
            water_signal = the_water_dictionary(:,t2w_i,flip_i);
            X = [water_signal sc_fat.signal];
            fs = X\signal;
            idx = sum(fs<0);
            distances(t2w_i,flip_i) = sum((signal-X*fs).^2)+idx*1e6; % Penalize negative fractions
            fractions(t2w_i,flip_i) = fs(1)/sum(fs);
        end
    end
    % end fit
    
    [COST(parind),IX] = min(distances(:));
    [xi,zi] = ind2sub(size(distances),IX);
    
    flip_angle_map(parind) = flip_angle_grid(zi);
    
    % refine the T2water while keeping all the other parameters constant
    distances = nan(length(T2w_finer_grid),1);
    fractions = nan(length(T2w_finer_grid),1);
    
    for t2w_i=1:length(T2w_finer_grid)
        water_signal = the_finer_water_dictionary(:,t2w_i,zi);
        X = [water_signal sc_fat.signal];
        fs = X\signal;
        idx = sum(fs<0);
        distances(t2w_i) = sum((signal-X*fs).^2)+idx*1e6;
        fractions(t2w_i) = fs(1)/sum(fs);
    end
    [COST(parind),IX] = min(distances(:));
    
    X=[the_finer_water_dictionary(:,IX,zi) sc_fat.signal];
    fs = lsqnonneg((X),(signal));
    
    T2w_map(parind) = T2w_finer_grid(IX);
    fw_map(parind) = fs(1)/sum(fs);
    ff_map(parind) = fs(2)/sum(fs);
    f_map(parind) = fractions(IX);
    % end refined T2water estimate
end
toc;

epg_fit.ff = zeros(size(T2SeriesData.img(:,:,:,1)));
epg_fit.fw = zeros(size(T2SeriesData.img(:,:,:,1)));
epg_fit.flip_angle_map = zeros(size(T2SeriesData.img(:,:,:,1)));
epg_fit.f_map = zeros(size(T2SeriesData.img(:,:,:,1)));
epg_fit.T2w_map = zeros(size(T2SeriesData.img(:,:,:,1)));
epg_fit.T2f_map = zeros(size(T2SeriesData.img(:,:,:,1)));
epg_fit.COST = zeros(size(T2SeriesData.img(:,:,:,1)));

epg_fit.ff(Points2Fit) = ff_map;
epg_fit.fw(Points2Fit) = fw_map;
epg_fit.flip_angle_map(Points2Fit) = flip_angle_map;
epg_fit.f_map(Points2Fit) = f_map;
epg_fit.T2w_map(Points2Fit) = T2w_map;
epg_fit.T2f_map(Points2Fit) = T2f_map;
epg_fit.COST(Points2Fit) = COST;

f_idx = strfind(T2MultiEchoNifti,'nii');
save([T2MultiEchoNifti(1:f_idx-1) '_T2EPG_BiExp'],'epg_fit');
end

% This function uses a late echo, where only fat is bright to obtain a
% subcutaneous fat segmentation
function Mask = GetSubcutaneousFat(T2SeriesData,echo)
    Mask = zeros(size(T2SeriesData.img(:,:,:,1)));
    eroder = strel('disk',3);
    for ij=1:size(Mask,3)
       Slice = squeeze(T2SeriesData.img(:,:,ij,echo)); 
       Slice = Slice/prctile(Slice(:),99);
       gt = prctile(Slice(:),90);%graythresh(Slice(:));
       Slice = Slice > gt;
       Mask(:,:,ij) = imerode(Slice,eroder);       
    end
    % Now find the largest connected component (subcutaneous fat)
    [L,N] = bwlabeln(Mask,26);
    NVox = zeros(N,1);
    for ij=1:N
       NVox(ij) = length(find(L==ij)); 
    end
    [~,IX] = max(NVox);
    for ij=1:N
        if(ij == IX)
            continue
        end
        Mask(L==ij) = 0;
    end
end

% EPG Matrix implementation of the Extended Phase Graph (EPG) algortihm.
%   H = EPG(N,TAU,R1,R2,ALPHA) calculate the echo amplitudes for a CPMG
%   sequence with N echoes and echo spacing of (2*TAU). The parameters are
%   R1=1/T1 (scalar) and R2=1./T2 (scalar or vector) and ALPHA the flip 
%   angle. The matrix H contains one column for each element in R2
%
% Author: Kelvin Layton (klayton@unimelb.edu.au)
% Date: June 2012
%
function [ H ] = epg( n, tau, R1, R2vec, alpha )


nRates = length(R2vec);
tau=tau/2;

H = zeros(n,nRates);

% RF mixing matrix
%
T0 = [cos(alpha/2)^2, sin(alpha/2)^2, sin(alpha); ...
    sin(alpha/2)^2, cos(alpha/2)^2, -sin(alpha); ...
    -0.5*sin(alpha), 0.5*sin(alpha), cos(alpha)];

TArray = cell(1,n);
TArray(:) = {sparse(T0)};
T = blkdiag(1,TArray{:});

% Selection matrix to move all traverse states up one coherence level
%
S = sparse(3*n+1,3*n+1);
S(1,3)=1;
S(2,1)=1;
S(3,6)=1;
S(4,4)=1;
for o=2:n
    offset1=( (o-1) - 1)*3 + 2;
    offset2=( (o+1) - 1)*3 + 3;
    if offset1<=(3*n+1)
    S(3*o-1,offset1)=1;  % F_k <- F_{k-1}
    end
    if offset2<=(3*n+1)
    S(3*o,offset2)=1;  % F_-k <- F_{-k-1}
    end
    S(3*o+1,3*o+1)=1;  % Z_order
end
    
for iRate=1:nRates

    % Relaxation matrix
    R2=R2vec(iRate);
    R0 = diag([exp(-tau*R2),exp(-tau*R2),exp(-tau*R1)]);

    RArray = cell(1,n);
    RArray(:) = {sparse(R0)};
    R = blkdiag(exp(-tau*R2),RArray{:});

    % Precession and relaxation matrix
    P = (R*S);

    % Matrix representing the inter-echo duration
    E = (P*T*P);
    
    % Recursively apply matrix to get echo amplitudes
    %
    x = zeros(size(R,1),1);
    x(1)=1;
    for iEcho=1:n
        x=E*x;
        H(iEcho,iRate) = x(1);
    end

end

end