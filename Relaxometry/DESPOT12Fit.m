%%%$ Included in MRIToolkit (https://github.com/delucaal/MRIToolkit) %%%%%% Alberto De Luca - alberto@isi.uu.nl $%%%%%% Distributed under the terms of LGPLv3  %%%
%%% Distributed under the terms of LGPLv3  %%%
% A. De Luca - UMC Utrecht - alberto@isi.uu.nl
% This script performs a DESPOT1/DESPOT2 fit from a balanced steady-state
% sequence (bSSFP) and a spoiled gradient echo acquisition (SPGR).

bSSFP = []; 
% prefix of the bSSFP sequence, assuming its name to be prefix_ph_01.nii
bSSFP_prefix = 'AD09_WIP_bSSFP_MF_20190603163703';

SPGR = [];
% prefix of the SPGR sequence, assuming its name to be prefix_ph_01.nii
SPGR_prefix = 'AD09_WIP_SPGR_MF_20190603163703';

FA_bSSFP = [22 18 14 10 6 2]; % The flip angles used for bSSFP
FA_SPGR = (2:2:18); % The flip angles used for SPGR

TR_bSSFP = 7.64; % The repetition time of the bSSFP
TR_SPGR = 8.85; % The repetition time of the SPGR

for ij=1:length(FA_bSSFP)
    bvol = load_untouch_nii([bSSFP_prefix '_ph_' sprintf('%02d',ij) '.nii']);
    bvol.img = single(bvol.img)*bvol.hdr.dime.scl_slope;%+bvol.hdr.dime.scl_inter;
    
    if(ij == 1)
        bSSFP = bvol.img;
    else
        bSSFP(:,:,:,end+1) = bvol.img;
    end
end

for ij=1:length(FA_SPGR)
    svol = load_untouch_nii([SPGR_prefix '_ph_' sprintf('%02d',ij) '.nii']);
    svol.img = single(svol.img)*svol.hdr.dime.scl_slope;%+svol.hdr.dime.scl_inter;
    
    if(ij == 1)
        SPGR = svol.img;
    else
        SPGR(:,:,:,end+1) = svol.img;
    end
end

% perform the DESPOT1 fit
T1_DESPOT1 = zeros(size(SPGR(:,:,:,1)));
for x=1:size(SPGR,1)
    for y=1:size(SPGR,2)
        for z=1:size(SPGR,3)
            S_SPGR = squeeze(SPGR(x,y,z,:))';
            if(any(S_SPGR == 0))
                continue
            end
            X = [ones(size(S_SPGR))' (S_SPGR./sind(FA_SPGR))'];
            p = X\(S_SPGR./tand(FA_SPGR))';
            slope = p(2);
            T1_DESPOT1(x,y,z) = abs(TR_SPGR/log(slope));
        end
    end
end

% perform the DESPOT2 fit
T2_DESPOT2 = zeros(size(bSSFP(:,:,:,1)));
for x=1:size(bSSFP,1)
    for y=1:size(bSSFP,2)
        for z=1:size(bSSFP,3)
            S_bSSFP = squeeze(bSSFP(x,y,z,:))';
            if(any(S_bSSFP == 0))
                continue
            end
            E1 = exp(-TR_bSSFP/T1_DESPOT1(x,y,z));
            X = [ones(size(S_bSSFP')) (S_bSSFP./sind(FA_bSSFP))']; 
            p = X\((S_bSSFP./tand(FA_bSSFP))');
            slope = p(2);
            T2_DESPOT2(x,y,z) = abs(TR_bSSFP/log(abs((slope-E1)/(slope*E1-1))));
        end
    end
end